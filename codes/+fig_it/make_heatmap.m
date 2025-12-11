function make_heatmap(p1,p2,f,p1_name,p2_name,plotOptions)

    valid = isfinite(p1) & isfinite(p2) & isfinite(f);
    p1 = p1(valid); p2 = p2(valid); f = f(valid);
    
    f_peak = zeros(size(f));
    p1_list = unique(p1);
    p2_list = unique(p2);
    for i = 1:numel(p1_list)
        for j = 1:numel(p2_list)
            idx = p1 == p1_list(i) & p2 == p2_list(j);
            f_min = min(f(idx)); f_max = max(f(idx));
            if abs(f_min) > abs(f_max)
                f_peak(idx) = f_min;
            else
                f_peak(idx) = f_max;
            end
           
        end
    end
    ng = 60; 
    xg = linspace(min(p1),max(p1),ng);
    yg = linspace(min(p2),max(p2),ng);
    [Xg,Yg] = meshgrid(xg,yg);
    if plotOptions.doPeak
        Fint = scatteredInterpolant(p1,p2,f_peak,'linear','linear'); % 'natural' is smooth
    else
        Fint = scatteredInterpolant(p1,p2,f,'linear','linear'); % 'natural' is smooth
    end
    Z = Fint(Xg,Yg);

    %nexttile; hold on
   
    contourf(Xg,Yg,Z,15,'LineStyle','none'); colorbar; hold on
   


    if isfield(plotOptions, 'f_obs')
        contour(Xg, Yg, Z, [plotOptions.f_obs plotOptions.f_obs], 'LineColor','k', 'LineWidth',2);
    end

    xlabel(p1_name); ylabel(p2_name); 
    set(gca, 'xtick',unique(p1));
    set(gca, 'ytick',unique(p2));
    set(gca,'fontsize',14)

end