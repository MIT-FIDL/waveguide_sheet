% 3D surface plot function
function [] = plot_curve_surf(app)
    plot_path = @(n_hplot,pts,color) plot_3D(app,n_hplot,pts,color);
    pos = @(endpt,t) [t*(endpt(3)-endpt(1))+endpt(1) t*(endpt(4)-endpt(2))+endpt(2)];

    evalc('colors = slanCM(''twilight'');');
    color_wg = @(dB) colors( min([max([int32(1+255*(dB+4.5)./9),1]) 255]), :);
    extention = 1+2/14;
    u = linspace(-extention,extention,app.N_mesh)./2;
    v = linspace(-extention,extention,app.N_mesh)./2;
    t = u'+0.5;
    n_u = [1 0 0]';
    n_v = [0 1 0]';
    R_mat = [n_u,n_v,cross(n_u,n_v)];
    P_mat = -app.P_uv(app.C,0,0);

    n_hplot = 1;
    for i = [-extention extention]/2      % plot perimeter
        pts = P_mat + R_mat*app.P_uv(app.C,u,i*ones(1,app.N_mesh));
        plot_path(n_hplot,pts,[]);n_hplot = n_hplot+1;
        pts = P_mat + R_mat*app.P_uv(app.C,i*ones(1,app.N_mesh),v);
        plot_path(n_hplot,pts,[]);n_hplot = n_hplot+1;
    end

    N_i = 1;
    for i = 1:app.N_wg      % plot waveguide path
        pts_uv = pos(app.wg_loc(:,i),t);
        pts = P_mat + R_mat*app.P_uv(app.C,pts_uv(:,1)',pts_uv(:,2)');
        plot_path(n_hplot,pts,color_wg(app.smoothedData(1,N_i)) );n_hplot = n_hplot+1;  
        pts_uv = pos(app.wg_loc(:,i+app.N_wg),t);
        pts = P_mat + R_mat*app.P_uv(app.C,pts_uv(:,1)',pts_uv(:,2)');
        plot_path(n_hplot,pts,color_wg(app.smoothedData(2,N_i)) );n_hplot = n_hplot+1;
        N_i = N_i+1; 
    end

    app.UIAxes.XLim = app.window_width;
    app.UIAxes.YLim = app.window_width;
    app.UIAxes.ZLim = app.window_width;
end
    
function [] = plot_3D(app,n_hplot,pts,colors)
    if ~isempty(colors)
        set(app.hPlot(n_hplot), 'XData', pts(1,:), 'YData', pts(2,:), 'ZData', pts(3,:), 'Color', colors);
    else
        set(app.hPlot(n_hplot), 'XData', pts(1,:), 'YData', pts(2,:), 'ZData', pts(3,:), 'Color', 'k');
    end
end
