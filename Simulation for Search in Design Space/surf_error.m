function [mean_error_avg,mean_error_std,max_error_avg,max_error_std] = surf_error(C_ref,N_wg,N_amp,tilt_ang,N_trial,noise_amp,plot_fig,s_n)

tilt_amp = 10*tand(tilt_ang/2);
mean_error_matrix = zeros(N_amp,N_trial);
max_error_matrix = zeros(N_amp,N_trial);
dB_per_cm_fit = @(kappa,n) kappa;

N_iteration = zeros(N_amp,1);
T = array2table([tilt_ang' sum(mean_error_matrix,2)./N_iteration sum(max_error_matrix,2)./N_iteration N_iteration],'VariableNames',{'Tilting Angle (deg)','Mean Norm Error','Max Norm Error','Iteration'}); 
disp(T);
warning ('off','all');

P_uv = @(A,u,v) A*[u; v; u.^2; u.*v; v.^2; u.^3; u.^2.*v; u.*v.^2; v.^3];
r_du = @(A,u,v) A*[ones(1,size(u,2));  zeros(1,size(u,2)); 2.*u; v; zeros(1,size(u,2)); 3.*u.^2; 2.*u.*v; v.^2;    zeros(1,size(u,2))];
r_dv = @(A,u,v) A*[zeros(1,size(u,2)); ones(1,size(u,2)); zeros(1,size(u,2)); u; 2.*v;  zeros(1,size(u,2)); u.^2; 2.*u.*v; 3.*v.^2];

for i = 1:N_amp
    wg_loc = get_wg_loc(N_wg,tilt_amp(i));
    color_wg = @(dB) [0.2 0.8 0.2];
    dB_ref = reverse_dB(C_ref,14,N_wg,wg_loc,dB_per_cm_fit);    
    rng(1816);
    for k = 1:N_trial
        dB = dB_ref + noise_amp*(rand(2,N_wg)-0.5);
        
        if plot_fig
            close all;hold on
            C = shape_optimize(N_wg,dB,wg_loc,plot_fig,color_wg);
            drawnow;
        else
            C = shape_optimize(N_wg,dB,wg_loc,plot_fig,color_wg);
        end

%         error_matrix(i,k) = get_norm_error(C_ref,C,r_du,r_dv);
%         error_matrix(i,k) = get_mad_dist_error(C_ref,C,P_uv);
        [mean_error_matrix(i,k),max_error_matrix(i,k)] = get_min_dist_error(C_ref,C,P_uv,r_du,r_dv);
        if isnan(mean_error_matrix(i,k))
            mean_error_matrix(i,k) = inf;
            warning('Failed to find minimum distance. The distance has been set to Inf.')
        end

        N_iteration(i) = N_iteration(i) + 1;
        T = array2table([tilt_ang' sum(mean_error_matrix,2)./N_iteration sum(max_error_matrix,2)./N_iteration N_iteration],'VariableNames',{'Tilting Angle (deg)','Mean Norm Error','Max Norm Error','Iteration'}); 
        clc;disp(T);
    end

    for k = 1:N_trial
        if mean_error_matrix(i,k) > 1e+10
            mean_error_matrix(i,k) = mean(mean_error_matrix(i,:)<=1e+10);
        end
    end
    
end

mean_error_avg = mean(mean_error_matrix,2);
mean_error_std = std(mean_error_matrix,1,2);
max_error_avg = mean(max_error_matrix,2);
max_error_std = std(max_error_matrix,1,2);
end

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

function norm_error = get_norm_error(C_ref,C,r_du,r_dv)
    norm_error = 0;
    N_grid_point = 10;
    uu = linspace(-0.5,0.5,N_grid_point); vv = uu;
    for i = 1:N_grid_point
        for j = 1:N_grid_point
            n1 = normalize( cross(r_du(C_ref,uu(i),vv(j)),r_dv(C_ref,uu(i),vv(j))), 'norm' );
            n2 = normalize( cross(r_du(C,uu(i),vv(j)),r_dv(C,uu(i),vv(j))),         'norm' );
            norm_error = norm_error + acos(dot(n1,n2));
        end
    end
    norm_error = norm_error/N_grid_point^2;
end

function dist_error = get_mad_dist_error(C_ref,C,P_uv)
    dist_error = 0;
    N_grid_point = 10;
    uu = linspace(-0.5,0.5,N_grid_point); vv = uu;

    R_ref = get_rot_matrix(C_ref(:,1)./norm(C_ref(:,1)), C_ref(:,2)./norm(C_ref(:,2)));
    R = get_rot_matrix(C(:,1)./norm(C(:,1)), C(:,2)./norm(C(:,2)));

    for i = 1:N_grid_point
        for j = 1:N_grid_point
            n1 = R_ref*P_uv(C_ref,uu(i),vv(j));
            n2 = R*P_uv(C,uu(i),vv(j));
            dist_error = dist_error + norm(n2-n1);
        end
    end
    dist_error = dist_error/N_grid_point^2;
end

function [mean_dist,max_dist] = get_min_dist_error(C_ref,C,P_uv,r_du,r_dv)
    mean_dist = 0;max_dist = 0;
    N_grid_point = 10;
    uu = linspace(-0.5,0.5,N_grid_point); vv = uu;

    R_ref = get_rot_matrix(C_ref(:,1)./norm(C_ref(:,1)), C_ref(:,2)./norm(C_ref(:,2)));
    R = get_rot_matrix(C(:,1)./norm(C(:,1)), C(:,2)./norm(C(:,2)));
    
    low = [-0.5;-0.5];up = [0.5;0.5];
    options.m = 5;         % Max L-BFGS pairs
    options.tol = 1e-6;    % Tolerance for convergence
    options.display = false; % Display iterations
    options.xhist = false;  % No need to store history in this case
%     distance_to_Sref = @(uv_on_ref,u1,v1) norm(R_ref*P_uv(C_ref,uv_on_ref(1),uv_on_ref(2)) - R*P_uv(C,u1,v1));

    for i = 1:N_grid_point
        for j = 1:N_grid_point
%             x0 = [uu(i);vv(j)];
            x0 = [uu(i);vv(j)];
            [x_min, ~] = LBFGSB(@(uv_on_ref) distance_to_Sref(uv_on_ref, uu(i), vv(j),C_ref,C,R_ref,R,P_uv,r_du,r_dv), x0, low, up, options);
            min_dist = distance_to_Sref(x_min, uu(i), vv(j), C_ref,C,R_ref,R,P_uv,r_du,r_dv);
            mean_dist = mean_dist + min_dist;
            max_dist = max(max_dist,min_dist);
        end
    end
    mean_dist = mean_dist/N_grid_point^2;
end

% Distance function and gradient between a point on S1 and a point on S2
function [dist, grad] = distance_to_Sref(uv_on_ref, u1, v1, C_ref,C,R_ref,R,P_uv,r_du,r_dv)
    % Extract parameters u2, v2
    u_ref = uv_on_ref(1);
    v_ref = uv_on_ref(2);

    % Compute points on surfaces S1 and S2
    P1 = R*P_uv(C,u1, v1);  % Point on S1
    Pref = R_ref*P_uv(C_ref,u_ref, v_ref);  % Point on S2

    % Compute distance (objective function)
    diff = P1 - Pref;
    dist = norm(diff);  % Euclidean distance

    % Compute gradient (partial derivatives of distance with respect to u2, v2)
    % Gradient of Euclidean distance:
    if dist > 0  % To avoid division by zero
        dP2_du2 = normalize(r_du(C_ref,u_ref, v_ref),'norm');  % Partial derivative of S2 with respect to u2
        dP2_dv2 = normalize(r_dv(C_ref,u_ref, v_ref),'norm');  % Partial derivative of S2 with respect to v2

        % Use the chain rule to compute the gradient
        grad_u2 = -diff' * dP2_du2 / dist;
        grad_v2 = -diff' * dP2_dv2 / dist;

        grad = [grad_u2; grad_v2];  % Gradient vector
    else
        grad = [0; 0];  % If distance is zero, gradient is zero
    end
end

function R = get_rot_matrix(n_u,n_v)
    R = [n_u n_v cross(n_u,n_v)]';
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

