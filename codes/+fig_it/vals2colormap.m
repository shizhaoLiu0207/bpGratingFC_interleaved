function rgb_colors = vals2colormap(vals, mapName, crange)
% Maps a vector of values to RGB colors from a specified colormap.
%
% Inputs:
% vals     = A vector of numeric values.
% mapName  = A string for the colormap name (e.g., 'jet', 'parula').
% crange   = [min_val, max_val] for color mapping (defaults to min/max of vals).
%
% Output:
% rgb_colors = Nx3 matrix of RGB colors [0, 1 range].
if nargin < 2
    mapName = 'turbo';
end

if nargin < 3
    crange = [min(vals), max(vals)]; % Default to full range of data
end

% Get the colormap matrix (default length is 256 colors)
cmap = feval(mapName, 256);
n_colors = size(cmap, 1);

% Scale the values to the colormap index range [1, n_colors]
% Handle potential edge cases for min/max values
scaled_vals = (vals - crange(1)) / (crange(2) - crange(1));
indices = round(scaled_vals * (n_colors - 1) + 1);

% Clamp indices to ensure they are within valid range [1, n_colors]
indices(indices < 1) = 1;
indices(indices > n_colors) = n_colors;

% Get the RGB colors using the calculated indices
rgb_colors = cmap(indices, :);

end