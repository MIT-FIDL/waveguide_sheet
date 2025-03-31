%% Initializing parameters
clear all;clc
N_amp = 20;         % Number of waveguide angle, higher N_amp corresponds to higher angle resolution
N_trial = 25;       % Number of trials per instance
noise_amp = 0.2;    % Percentage of noise (total range, equal to plus or minus half of noise_amp)
tilt_ang = logspace(0,log10(19),N_amp)-1;
N_wg_arr = 3:20;
N_wg_size = size(N_wg_arr,2);
P_uv = @(A,u,v) A*[u; v; u.^2; u.*v; v.^2; u.^3; u.^2.*v; u.*v.^2; v.^3];
r_du = @(A,u,v) A*[ones(1,size(u,2));  zeros(1,size(u,2)); 2.*u; v; zeros(1,size(u,2)); 3.*u.^2; 2.*u.*v; v.^2;    zeros(1,size(u,2))];
r_dv = @(A,u,v) A*[zeros(1,size(u,2)); ones(1,size(u,2)); zeros(1,size(u,2)); u; 2.*v;  zeros(1,size(u,2)); u.^2; 2.*u.*v; 3.*v.^2];

%% Shape S1
% C_ref = reference C from optimizaion
C_ref = [13.3597   -3.2348    0.2151    0.6482    0.5584   -1.3715   -3.9529   -3.4152   -2.9481;...
         2.4510    9.5794    4.8402    7.6542    4.1323   -2.1511   -3.5997   -4.4498   -2.0871;...
        -3.8927   -9.8860    3.6356    7.2145    3.6585    2.6019    5.2012    4.5370    2.0604];
plot_ref_surface(C_ref,get_wg_loc(10,0),[-0.2 -1 0.3]); % plot reference shape

mean_error_avg = zeros(N_amp,N_wg_size);
mean_error_std = zeros(N_amp,N_wg_size);
max_error_avg = zeros(N_amp,N_wg_size);
max_error_std = zeros(N_amp,N_wg_size);
for k = 1:N_wg_size % loop through design space for distance error of shape reconstruction
    [mean_error_avg(:,k),mean_error_std(:,k),max_error_avg(:,k),max_error_std(:,k)] = surf_error(C_ref,N_wg_arr(k),N_amp,tilt_ang,N_trial,noise_amp,true);
end
save('surf_error_analysis/s_error_1_noise_10.mat','-mat')

%% Shape S2
C_ref = [13.2479   -3.9161   -0.0671   -0.4797    0.0027   -0.6746    0.2338   -0.8924    0.6529;...
         4.1086   13.1700    0.0804    0.1256   -0.1509   -0.1431    1.0075    0.0442   -0.0985;...
         1.8199   -0.2430    2.7574   -4.0286    1.7748    0.5386    7.2646   -6.0659   -2.8762];
plot_ref_surface(C_ref,get_wg_loc(7,0),[-0.3 -1 0.4])

mean_error_avg = zeros(N_amp,N_wg_size);
mean_error_std = zeros(N_amp,N_wg_size);
max_error_avg = zeros(N_amp,N_wg_size);
max_error_std = zeros(N_amp,N_wg_size);
for k = 1:N_wg_size
    [mean_error_avg(:,k),mean_error_std(:,k),max_error_avg(:,k),max_error_std(:,k)] = surf_error(C_ref,N_wg_arr(k),N_amp,tilt_ang,N_trial,noise_amp,true,P_uv,r_du,r_dv);
end
save('surf_error_analysis/s_error_2_noise_10.mat','-mat')

%% Shape S3
C_ref = [13.3554   -1.6626   -0.1581    0.2039   -0.0162   -3.1756   7.0363    2.6143   -2.7252;...
         1.9734   13.9643   -0.0830   -0.1172    0.0373   -1.6590    1.8341   -0.8316   -1.2255;...
        -5.0262   -0.1681   -0.0190    0.1992   -0.1162   -4.8525   13.9862   12.0213   -1.5681];
plot_ref_surface(C_ref,get_wg_loc(7,0),[-0.3 -1 0.3])

mean_error_avg = zeros(N_amp,N_wg_size);
mean_error_std = zeros(N_amp,N_wg_size);
max_error_avg = zeros(N_amp,N_wg_size);
max_error_std = zeros(N_amp,N_wg_size);
for k = 1:N_wg_size
    [mean_error_avg(:,k),mean_error_std(:,k),max_error_avg(:,k),max_error_std(:,k)] = surf_error(C_ref,N_wg_arr(k),N_amp,tilt_ang,N_trial,noise_amp,true,P_uv,r_du,r_dv);
end
save('surf_error_analysis/s_error_3_noise_10.mat','-mat')

%% Shape S4
C_ref = [12.9166    1.7268    0.0621    0.0802    1.5385    0.0391   -0.3271   -4.5355    0.0052;...
         -2.0097   14.0134    0.0172    1.5199    0.1450   -0.0375   -0.7067   -0.1571   -1.4376;...
         -5.3991   -1.0658    0.0547   -0.3204    4.0722    0.0656   -0.4825   -9.3968    0.0446];
plot_ref_surface(C_ref,get_wg_loc(7,0),[1 -0.3 0.4])

mean_error_avg = zeros(N_amp,N_wg_size);
mean_error_std = zeros(N_amp,N_wg_size);
max_error_avg = zeros(N_amp,N_wg_size);
max_error_std = zeros(N_amp,N_wg_size);
for k = 1:N_wg_size
    [mean_error_avg(:,k),mean_error_std(:,k),max_error_avg(:,k),max_error_std(:,k)] = surf_error(C_ref,N_wg_arr(k),N_amp,tilt_ang,N_trial,noise_amp,true,P_uv,r_du,r_dv);
end
save('surf_error_analysis/s_error_4_noise_10.mat','-mat')

%% Shape S5
C_ref = [14.0816   -0.1500    0.0770    0.6202   -0.1630   -1.3087    0.6800   -1.4288    0.2088;...
        -0.1211   14.0218   -0.1386    0.0567    0.0248    1.0870   -1.2155    0.4827   -0.0923;...
         0.5490   -0.6439   -0.6311    0.5147    0.1561    0.3389   17.3616   -3.9429    0.6541];
plot_ref_surface(C_ref,get_wg_loc(7,0),[-0.45 -1 0.45])

mean_error_avg = zeros(N_amp,N_wg_size);
mean_error_std = zeros(N_amp,N_wg_size);
max_error_avg = zeros(N_amp,N_wg_size);
max_error_std = zeros(N_amp,N_wg_size);
for k = 1:N_wg_size
    [mean_error_avg(:,k),mean_error_std(:,k),max_error_avg(:,k),max_error_std(:,k)] = surf_error(C_ref,N_wg_arr(k),N_amp,tilt_ang,N_trial,noise_amp,true,P_uv,r_du,r_dv);
end
save('surf_error_analysis/s_error_5_noise_10.mat','-mat')

%% Shape S6
C_ref = [13.0557    0.4111   -6.3093   -0.0054   -0.0169  -14.3752   -0.0195   -0.0209   -0.1778;...
         0.4621   14.0061    2.0509    0.0019    0.0018   -0.0699   -0.0083   -0.0057    0.0745
         6.1857   -1.9351   13.5327    0.0084   -0.0269   -7.0590   -0.0070   -0.0202    0.0240];
plot_ref_surface(C_ref,get_wg_loc(7,0),[-0.4 -1 0.35])

mean_error_avg = zeros(N_amp,N_wg_size);
mean_error_std = zeros(N_amp,N_wg_size);
max_error_avg = zeros(N_amp,N_wg_size);
max_error_std = zeros(N_amp,N_wg_size);
for k = 1:N_wg_size
    [mean_error_avg(:,k),mean_error_std(:,k),max_error_avg(:,k),max_error_std(:,k)] = surf_error(C_ref,N_wg_arr(k),N_amp,tilt_ang,N_trial,noise_amp,true,P_uv,r_du,r_dv);
end
save('surf_error_analysis/s_error_6_noise_10.mat','-mat')

%% Heatmap ploting
load('surf_error_analysis/s_error_1_noise_10.mat','mean_error_avg','tilt_ang','N_wg_arr')
[X, Y] = meshgrid([tilt_ang tilt_ang(end)+3.35],[N_wg_arr N_wg_arr(end)+1]);
Y = Y-0.5;X = X-0.5;
D1 = normalize_D(mean_error_avg);
load('surf_error_analysis/s_error_2_noise_10.mat','mean_error_avg')
D2 = normalize_D(mean_error_avg);
load('surf_error_analysis/s_error_3_noise_10.mat','mean_error_avg')
D3 = normalize_D(mean_error_avg);
load('surf_error_analysis/s_error_4_noise_10.mat','mean_error_avg')
D4 = normalize_D(mean_error_avg);
load('surf_error_analysis/s_error_5_noise_10.mat','mean_error_avg')
D5 = normalize_D(mean_error_avg);
load('surf_error_analysis/s_error_6_noise_10.mat','mean_error_avg')
D6 = normalize_D(mean_error_avg);

D = D1+D2+D3+D5+D6;
D_max = max(max(D));
D_min = 0;
D = (D-D_min)./(D_max-D_min);
figure;hold on
pcolor(X,Y, D7);
shading flat;
pbaspect([1 1 1])
colormap(slanCM('Oranges'))
cb = colorbar;
clim([0 0.7879]);
yl = ylabel(cb,'Normalized Mean Distance to Target Shape','FontSize',12,'Rotation',270);

set(gca, 'FontSize', 12);
ylabel('# Waveguides per Axis');
xlabel('\theta (Deg)');
xlim([-0.5 20.85])
ylim([2.5 20.5])

%% Functions
function C = shape_optimize(N_wg,sensor_dB,wg_loc,plot_surf,color_wg)
L = 14;
N_coeff = 9;
C = zeros(3,N_coeff);C(1,1) = L-0.2; C(2,2) = L-0.2;

dB_per_cm_fit = @(kappa,n) kappa;
options = optimoptions('fmincon','Algorithm','sqp','EnableFeasibilityMode',true,'SubproblemAlgorithm', 'cg',...
    'MaxIterations',1e+5,'SpecifyObjectiveGradient',false,'SpecifyConstraintGradient',false,...
    'ConstraintTolerance',5,'FunctionTolerance',1e-5,'StepTolerance',1e-5,'Display','none','MaxFunctionEvaluations',1e+5,...
    'ObjectiveLimit',1,'OptimalityTolerance',1e-5);
P_uv = @(A,u,v) A*[u; v; u.^2; u.*v; v.^2; u.^3; u.^2.*v; u.*v.^2; v.^3];

C = shape_optimizer(sensor_dB,N_coeff,C,N_wg,wg_loc,L,dB_per_cm_fit,options,false,[1 0.3 0.5]);

if plot_surf
    plot_curve_surf(C,P_uv,N_wg,wg_loc,sensor_dB,[0.3 -1 0.5],color_wg)
end

end

function wg_loc = get_wg_loc(N_wg,tilt_amp)
wg_loc = zeros(4,2*N_wg);   % first 2 row = u1,v1, last 2 row = u2,v2, wg 1-7 and then 8-14
wg_loc(1,1:N_wg) = -0.5; wg_loc(3,1:N_wg) = 0.5;
wg_loc([2 4],1:N_wg) = [linspace(-0.45,0.45,N_wg);linspace(-0.45,0.45,N_wg)];
wg_loc(2,1+N_wg:end) = -0.5; wg_loc(4,1+N_wg:end) = 0.5;
wg_loc([1 3],1+N_wg:end) = [linspace(-0.45,0.45,N_wg);linspace(-0.45,0.45,N_wg)];

wg_loc(2,2:2:N_wg) = wg_loc(2,2:2:N_wg)-0.05*tilt_amp;
wg_loc(4,2:2:N_wg) = wg_loc(4,2:2:N_wg)+0.05*tilt_amp;
wg_loc(1,(2:2:N_wg)+N_wg) = wg_loc(1,(2:2:N_wg)+N_wg)+0.05*tilt_amp;
wg_loc(3,(2:2:N_wg)+N_wg) = wg_loc(3,(2:2:N_wg)+N_wg)-0.05*tilt_amp;
wg_loc(2,1:2:N_wg) = wg_loc(2,1:2:N_wg)+0.05*tilt_amp;
wg_loc(4,1:2:N_wg) = wg_loc(4,1:2:N_wg)-0.05*tilt_amp;
wg_loc(1,(1:2:N_wg)+N_wg) = wg_loc(1,(1:2:N_wg)+N_wg)-0.05*tilt_amp;
wg_loc(3,(1:2:N_wg)+N_wg) = wg_loc(3,(1:2:N_wg)+N_wg)+0.05*tilt_amp;
wg_loc(:,8:end) = flip(wg_loc(:,8:end),2);
end

function max_tilt_ang = plot_design_space(N_wg_arr,plot_flag)
w = 6;  % waveguide width
L_wg = 162;
max_tilt_ang = zeros(1,size(N_wg_arr,2));
syms theta
L = L_wg*cosd(theta) + w*sind(theta);
L_vert = N_wg_arr*(L_wg*sind(theta) + w*cosd(theta));

for i = 1:size(N_wg_arr,2)
    ang = double(solve(L == L_vert(i),theta));
    max_tilt_ang(i) = real(ang(ang>=0));
end
if plot_flag
    hold on;
    plot(max_tilt_ang,N_wg_arr,'k--','LineWidth',1.5)
end

end

function plot_ref_surface(C_mat,wg_loc,view_vec)
P_uv = @(A,u,v) A*[u; v; u.^2; u.*v; v.^2; u.^3; u.^2.*v; u.*v.^2; v.^3];
% P_uv = @(A,u,v) A*[u; v; -(10.*u.^2.*v - u.^2 - v.^2 + 30*u.^3.*v.^2 - v.^4 - 1.5.*u - v)./5];

evalc('colors = colormap(slanCM(''twilight''));');
color_wg = @(dB) [0.2 0.8 0.2];

plot_curve_surf(C_mat,P_uv,size(wg_loc,2)/2,wg_loc,ones(2,10),view_vec,color_wg)
xlim([-15 15])
ylim([-15 15])
zlim([-5 15])
end

function D = normalize_D(data_matrix)
size_data = size(data_matrix);
data_matrix = smoothdata2(data_matrix,'lowess',2.1);
data_matrix = [data_matrix zeros(size_data(1),1)];
data_matrix = [data_matrix;zeros(1,size_data(2)+1)];

D = data_matrix';
end