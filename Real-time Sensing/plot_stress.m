function mag = plot_stress(app)

pxyz = app.P_uv(app.C,app.stress_pts(1),app.stress_pts(2));
mag = 20*( max(max(app.smoothedData))-min(min(app.smoothedData)) );
% hold on;
% plot3(pxyz(1),pxyz(2),pxyz(3),'r.','MarkerSize',20*(max(max(sensor_dB))-min(val1)) )
if mag >= 10
    set(app.hPlot(end), 'XData', pxyz(1), 'YData', pxyz(2), 'ZData', pxyz(3), 'MarkerSize', mag);
else
    set(app.hPlot(end), 'XData', 100, 'YData', 100, 'ZData', 100, 'MarkerSize', 0.1);
end

end
