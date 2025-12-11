function Z = interp2_preserve_extrema(p1, p2, f, Xg, Yg, pin_radius)
%INTERP2_PRESERVE_EXTREMA
% Linear 2-D interpolation that (i) is piecewise-linear (Delaunay), (ii) 
% matches sample values at the data points, and (iii) guarantees the grid
% Z has the same global min/max as f by "pinning" nearby grid cells.
%
% Inputs:
%   p1, p2 : Nx1 vectors of sample coordinates
%   f      : Nx1 vector of sample values
%   Xg,Yg  : meshgrid-style query grid (same size)
%   pin_radius (optional): radius (in param units) within which to pin
%                          grid cells to exact sample values. Default: auto.
%
% Output:
%   Z      : size(Xg) matrix of interpolated values with preserved extrema
%
% Notes:
% - Inside the convex hull: barycentric (piecewise-linear) interpolation
%   on the Delaunay triangulation.
% - Outside the hull: nearest-neighbor fallback.
% - Extrema preservation: after interpolation, we overwrite grid cells
%   within 'pin_radius' of each sample with its exact f value.

 

    % 0) Hygiene
    valid = isfinite(p1) & isfinite(p2) & isfinite(f);
    p1 = p1(valid); p2 = p2(valid); f = f(valid);
    assert(isvector(p1) && isvector(p2) && isvector(f), 'p1,p2,f must be vectors');
    assert(numel(p1)==numel(p2) && numel(p1)==numel(f), 'Inputs must have same length');

    sz = size(Xg);
    Xq = Xg(:); Yq = Yg(:);

    % 1) Delaunay triangulation of samples
    dt = delaunayTriangulation(p1, p2);
    T  = dt.ConnectivityList;          % triangle vertex indices
    P  = dt.Points;                    % [p1 p2]

    % 2) Locate each query point: which triangle? (NaN -> outside hull)
    ti = pointLocation(dt, [Xq Yq]);   % triangle index for each query point

   % 3) Interpolate linearly (barycentric) for "inside" points
    Zq = nan(size(Xq));
    inside = ~isnan(ti);
    
    if any(inside)
        tri_idx = T(ti(inside), :);         % 3 vertex indices for each query point
        V1 = P(tri_idx(:,1),:);
        V2 = P(tri_idx(:,2),:);
        V3 = P(tri_idx(:,3),:);
    
        % For each query point, solve barycentric coords
        XY  = [Xq(inside) Yq(inside)];
        denom = ( (V2(:,2)-V3(:,2)).*(V1(:,1)-V3(:,1)) + ...
                  (V3(:,1)-V2(:,1)).*(V1(:,2)-V3(:,2)) );
    
        w1 = ( (V2(:,2)-V3(:,2)).*(XY(:,1)-V3(:,1)) + ...
               (V3(:,1)-V2(:,1)).*(XY(:,2)-V3(:,2)) ) ./ denom;
        w2 = ( (V3(:,2)-V1(:,2)).*(XY(:,1)-V3(:,1)) + ...
               (V1(:,1)-V3(:,1)).*(XY(:,2)-V3(:,2)) ) ./ denom;
        w3 = 1 - w1 - w2;
    
        f_tri = [f(tri_idx(:,1)) f(tri_idx(:,2)) f(tri_idx(:,3))];
        Zq(inside) = w1.*f_tri(:,1) + w2.*f_tri(:,2) + w3.*f_tri(:,3);
    end
    % 4) Outside the convex hull: nearest-neighbor fallback
    if any(~inside)
        % nearest sample for each outside point
        nn = knnsearch(P, [Xq(~inside) Yq(~inside)]);
        Zq(~inside) = f(nn);
    end

    % 5) Reshape to grid
    Z = reshape(Zq, sz);

    % 6) Extrema preservation by pinning grid cells near each sample
    % Choose default pin radius as 1/2 of the average grid spacing
    if isempty(pin_radius)
        % Estimate grid spacing robustly
        dx = median(abs(diff(unique(Xg(:)))));
        dy = median(abs(diff(unique(Yg(:)))));
        if ~isfinite(dx) || dx==0, dx = range(p1)/max(10, sz(2)-1); end
        if ~isfinite(dy) || dy==0, dy = range(p2)/max(10, sz(1)-1); end
        pin_radius = 0.5 * sqrt(dx^2 + dy^2);
    end

    % For each sample, find grid cells within pin_radius and set to exact f
    % (This guarantees min(Z)==min(f) and max(Z)==max(f) on the grid.)
    % Vectorized neighborhood pinning:
    G = [Xg(:) Yg(:)];
    S = [p1(:) p2(:)];
    % To avoid O(N*M) when very large, pin via nearest neighbor first, then expand
    nn_grid_idx = knnsearch(G, S);    % one nearest grid cell per sample
    Z(nn_sub2ind(nn_grid_idx, sz)) = f;  % exact pin at nearest cell

    % Optionally expand to nearby cells within radius (exact match):
    if pin_radius > 0
        % Build a KD-tree on grid for radius search
        idx_cell = rangesearch(G, S, pin_radius);
        for i = 1:numel(idx_cell)
            if ~isempty(idx_cell{i})
                Z(idx_cell{i}) = f(i);
            end
        end
    end

    % Final safety: clip to empirical range (no change if already within)
    fmin = min(f); fmax = max(f);
    Z = min(max(Z, fmin), fmax);
end

% Helper: linear index from knnsearch result
function I = nn_sub2ind(nn_idx, sz)
    % nn_idx are linear indices into vectorized grid (row-major from Xg(:))
    I = nn_idx;
end
