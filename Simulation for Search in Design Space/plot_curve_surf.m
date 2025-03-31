function [] = plot_curve_surf(C_mat,P_uv,N_wg,wg_loc,sensor_dB,view_angle,color_wg)
hold on
N_mesh = 30;
extention = 1+2/14;
u = linspace(-extention,extention,N_mesh)./2;
v = linspace(-extention,extention,N_mesh)./2;
t = u'+0.5;
% uv_vector = @(u,v) [u; v; u.^2; u.*v; v.^2];

plot_path = @(pts,linestyle,linewidth,color) plot_3D(pts,linestyle,linewidth,color);
pos = @(endpt,t) [t*(endpt(3)-endpt(1))+endpt(1) t*(endpt(4)-endpt(2))+endpt(2)];

% % n_u = -P_uv(C_mat,0.001,0)./norm(P_uv(C_mat,0.001,0));
% % n_v = -P_uv(C_mat,0,0.001)./norm(P_uv(C_mat,0,0.001));
% n_u = C_mat(:,1)./norm(C_mat(:,1));
% n_v = C_mat(:,2)./norm(C_mat(:,2));
% n_norm = cross(n_u,n_v);
% n_rot = normalize(cross(n_norm,[0 0 1]'),'norm');
% phi = acos(n_norm(3));
% K_cross = [0 -n_rot(3) n_rot(2);n_rot(3) 0 -n_rot(1);-n_rot(2) n_rot(1) 0];
% R_mat = eye(3)+sin(phi)*K_cross + (1-cos(phi))*K_cross^2;
% n_u_planar = normalize(n_u(1:2),'norm');
% theta = -sign(n_u_planar(2))*acos(n_u_planar(1));
% % theta = 0;
% R_mat = [cos(theta) -sin(theta) 0;sin(theta) cos(theta) 0;0 0 1]*R_mat;
% % R_mat = eye(3)+ [(1-cos(angle))*(n_rot(1)^2-1) -n_rot(3)*sin(angle)+(1-cos(angle))*n_rot(1)*n_rot(2) n_rot(2)*sin(angle)+(1-cos(angle))*n_rot(1)*n_rot(3);...
% %                  n_rot(3)*sin(angle)+(1-cos(angle))*n_rot(1)*n_rot(2) (1-cos(angle))*(n_rot(2)^2-1) -n_rot(1)*sin(angle)+(1-cos(angle))*n_rot(2)*n_rot(3);...
% %                  -n_rot(2)*sin(angle)+(1-cos(angle))*n_rot(1)*n_rot(3) n_rot(1)*sin(angle)+(1-cos(angle))*n_rot(2)*n_rot(3) (1-cos(angle))*(n_rot(3)^2-1)];
% 
% if isnan(n_rot(1))
%     R_mat = eye(3);
% end

n_u = C_mat(:,1)./norm(C_mat(:,1));
n_v = C_mat(:,2)./norm(C_mat(:,2));
R_mat = [n_u n_v cross(n_u,n_v)]';

% t_x = -12;
% t_y = 10;
% t_z = 0;
% % t_x = -12;
% % t_y = -2;
% % t_z = 1;
% % t_x = 6;
% % t_y = -16;
% % t_z = 4;
% K_cross = [0 0 0;0 0 -1;0 1 0];
% R_mat = (eye(3)+sind(t_x)*K_cross + (1-cosd(t_x))*K_cross^2) * R_mat;
% K_cross = [0 0 1;0 0 0;-1 0 0];
% R_mat = (eye(3)+sind(t_y)*K_cross + (1-cosd(t_y))*K_cross^2) * R_mat;
% K_cross = [0 -1 0;1 0 0;0 0 0];
% R_mat = (eye(3)+sind(t_z)*K_cross + (1-cosd(t_z))*K_cross^2) * R_mat;
T_mat = [0 0 0]'; % 0 0 1 for s3

N_mesh_f = 30;
extention = 1+2/14;
uu = linspace(-extention,extention,N_mesh_f+1)./2;
vv = linspace(-extention,extention,N_mesh_f+1)./2;
for i = 1:N_mesh_f
    for j = 1:N_mesh_f
        pts = meshgrid(T_mat,ones(4,1))' + R_mat*[P_uv(C_mat,uu(i),vv(j)) P_uv(C_mat,uu(i),vv(j+1)) P_uv(C_mat,uu(i+1),vv(j+1)) P_uv(C_mat,uu(i+1),vv(j))];
        fill3(pts(1,:),pts(2,:),pts(3,:),'w','EdgeColor','none','FaceAlpha',1,'FaceColor',1*[1 1 1])
    end
end

N_i = 1;
for i = 1:N_wg      % plot waveguide path
    pts_uv = pos(wg_loc(:,i),t);
    pts = meshgrid(T_mat,ones(N_mesh,1))' + R_mat*P_uv(C_mat,pts_uv(:,1)',pts_uv(:,2)');
    plot_path(pts,'-',1.5,color_wg(sensor_dB(1,N_i)) )    
    pts_uv = pos(wg_loc(:,i+N_wg),t);
    pts = meshgrid(T_mat,ones(N_mesh,1))' + R_mat*P_uv(C_mat,pts_uv(:,1)',pts_uv(:,2)');
    plot_path(pts,'-',1.5,color_wg(sensor_dB(2,N_i)) )    
    N_i = N_i+1;
end

for i = [-extention extention]./2      % plot perimeter
    pts = meshgrid(T_mat,ones(N_mesh,1))' + R_mat*P_uv(C_mat,u,i*ones(1,N_mesh));
    plot_path(pts,'k-',3,[])
    pts = meshgrid(T_mat,ones(N_mesh,1))' + R_mat*P_uv(C_mat,i*ones(1,N_mesh),v);
    plot_path(pts,'k-',3,[])    
end

view(view_angle);
grid on; axis equal
xlim([-10 10]);ylim([-10 10]);zlim([-10 10])
ax = gca;
ax.FontSize = 14; 
xticks(-12:4:12)
yticks(-12:4:12)
zticks(-8:4:8)

% xlabel('x (cm)', 'FontSize', 14)
% ylabel('y (cm)', 'FontSize', 14)
% zlabel('z (cm)', 'FontSize', 14)
end

function [] = plot_3D(pts,linestyle,linewidth,colors)
if ~isempty(colors)
    plot3(pts(1,:),pts(2,:),pts(3,:),linestyle,'LineWidth',linewidth,'Color',colors)
else
    plot3(pts(1,:),pts(2,:),pts(3,:),linestyle,'LineWidth',linewidth)
end
end
