function pts = stress_optimizer(sensor_dB,P_uv,C,wg_loc,plot_result)
signal = sensor_dB;
[val1,ind1] = min(signal,[],2);
signal(1,ind1(1)) = inf;
signal(2,ind1(2)) = inf;
[val2,ind2] = min(signal,[],2);

w_x = val1(1)^2 / (val1(1)^2 + val2(1)^2);
w_y = val1(2)^2 / (val1(2)^2 + val2(2)^2);
min_loc = wg_loc(:,[ind1(1) ind2(1)]);
min_loc = w_x*min_loc(:,1) + (1-w_x)*min_loc(:,2);
uv_x = [min_loc([1 3])'; min_loc([2 4])'];
min_loc = wg_loc(:,7+[ind1(2) ind2(2)]);
min_loc = w_y*min_loc(:,1) + (1-w_y)*min_loc(:,2);
uv_y = [min_loc([1 3])'; min_loc([2 4])'];

a_x = (uv_x(2,2)-uv_x(2,1))/(uv_x(1,2)-uv_x(1,1));
b_x = uv_x(2,1) - (uv_x(2,2)-uv_x(2,1))/(uv_x(1,2)-uv_x(1,1))*uv_x(1,1);
a_y = (uv_y(2,2)-uv_y(2,1))/(uv_y(1,2)-uv_y(1,1));
b_y = uv_y(2,1) - (uv_y(2,2)-uv_y(2,1))/(uv_y(1,2)-uv_y(1,1))*uv_y(1,1);

u_str = (b_y-b_x)/(a_x-a_y);
pts = [u_str; a_x*u_str+b_x];
% pts = val1.^2./(val1.^2+val2.^2).*wg_loc(ind1)' + val2.^2./(val1.^2+val2.^2).*wg_loc(ind2)';


if plot_result
    evalc('colors = colormap(slanCM(''twilight''));');
    color_wg = @(dB) colors( min([max([int32(1+255*(dB+4.5)./9),1]) 255]), :);
    pxyz = P_uv(C,pts(1),pts(2));
    hold on;
    plot_curve_surf(C,P_uv,7,wg_loc,sensor_dB,[1 0.5 0.5],color_wg)
    plot3(pxyz(1),pxyz(2),pxyz(3),'r.','MarkerSize',20*(max(max(sensor_dB))-min(val1)) )
    cb = colorbar;caxis([-4.5 4.5]);yl = ylabel(cb,'Power (dB)','FontSize',14,'Rotation',270);yl.Position(1) = yl.Position(1)+1;
    cb.Position(1) = cb.Position(1) +0.03;
    hold off
end
end
