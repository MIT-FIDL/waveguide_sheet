function [C_mat,fval] = shape_optimizer(sensor_dB,N_coeff,C0,N_wg,wg_loc,L,dB_per_cm_fit,options,plot_result,view_vec)
if N_coeff ~= 9
    error('No. Coefficient Mismatch')
end

constraint_eval_loc = (1:4:13)./14-0.5;
% constraint_eval_loc = linspace(-0.45,0.45,10);

C0 = C0(:);

%% Specify P,dP/du,dP/dv,ddP/duu,ddP/dvv

% r_du = @(A,u,v) A*[ones(1,size(u,2));  zeros(1,size(u,2)); 2.*u; v; zeros(1,size(u,2)); 3.*u.^2; 2.*u.*v; v.^2;    zeros(1,size(u,2));...
%                    2.*u.*v.^2; 3.*u.^2.*v; v.^3; 3.*u.^2.*v.^2; 2.*u.*v.^3; 3.*u.^2.*v.^3];
% r_dv = @(A,u,v) A*[zeros(1,size(u,2)); ones(1,size(u,2)); zeros(1,size(u,2)); u; 2.*v;  zeros(1,size(u,2)); u.^2; 2.*u.*v; 3.*v.^2;...
%                    2.*u.^2.*v; u.^3; 3.*u.*v.^2; 2.*u.^3.*v; 3.*u.^2.*v.^2; 3.*u.^3.*v.^2];
% 
% r_ddu = @(A,u,v) A(:,[3 6 7 10:end])*[2*ones(1,size(u,2)); 6.*u; 2.*v; ...
%                    2.*v.^2; 6.*u.*v; zeros(1,size(u,2)); 6.*u.*v.^2; 2.*v.^3; 6.*u.*v.^3];
% r_dudv= @(A,u,v) A(:,[4 7 8 10:end])*[ones(1,size(u,2)); 2.*u; 2.*v; ...
%                    4.*u.*v; 3.*u.^2; 3.*v.^2; 6.*u.^2.*v; 6.*u.*v.^2; 9.*u.^2.*v.^2];
% r_ddv = @(A,u,v) A(:,[5 8 9 10:end])*[2*ones(1,size(u,2)); 2.*u; 6.*v; ...
%                    2.*u.^2; zeros(1,size(u,2)); 6.*u.*v; 2.*u.^3; 6.*u.^2.*v; 6.*u.^3.*v];

% r_du = @(A,u,v) A*[ones(1,size(u,2));  zeros(1,size(u,2)); 2.*u; v; zeros(1,size(u,2)); 3.*u.^2; 2.*u.*v; v.^2;    zeros(1,size(u,2));...
%                    3.*u.^2.*v.^2; 2.*u.*v.^3; 3.*u.^2.*v.^3];
% r_dv = @(A,u,v) A*[zeros(1,size(u,2)); ones(1,size(u,2)); zeros(1,size(u,2)); u; 2.*v;  zeros(1,size(u,2)); u.^2; 2.*u.*v; 3.*v.^2;...
%                    2.*u.^3.*v; 3.*u.^2.*v.^2; 3.*u.^3.*v.^2];
% 
% r_ddu = @(A,u,v) A(:,[3 6 7 10:end])*[2*ones(1,size(u,2)); 6.*u; 2.*v; ...
%                    6.*u.*v.^2; 2.*v.^3; 6.*u.*v.^3];
% r_dudv= @(A,u,v) A(:,[4 7 8 10:end])*[ones(1,size(u,2)); 2.*u; 2.*v; ...
%                    6.*u.^2.*v; 6.*u.*v.^2; 9.*u.^2.*v.^2];
% r_ddv = @(A,u,v) A(:,[5 8 9 10:end])*[2*ones(1,size(u,2)); 2.*u; 6.*v; ...
%                    2.*u.^3; 6.*u.^2.*v; 6.*u.^3.*v];

r_du = @(A,u,v) A*[ones(1,size(u,2));  zeros(1,size(u,2)); 2.*u; v; zeros(1,size(u,2)); 3.*u.^2; 2.*u.*v; v.^2;    zeros(1,size(u,2))];
r_dv = @(A,u,v) A*[zeros(1,size(u,2)); ones(1,size(u,2)); zeros(1,size(u,2)); u; 2.*v;  zeros(1,size(u,2)); u.^2; 2.*u.*v; 3.*v.^2];
r_ddu = @(A,u,v) A(:,[3 6 7])*[2*ones(1,size(u,2)); 6.*u; 2.*v];
r_dudv= @(A,u,v) A(:,[4 7 8])*[ones(1,size(u,2)); 2.*u; 2.*v];
r_ddv = @(A,u,v) A(:,[5 8 9])*[2*ones(1,size(u,2)); 2.*u; 6.*v];

%% Run fmincon
% lb = -Inf*ones(3*N_coeff,1);
% ub = Inf*ones(3*N_coeff,1);
% lb(end-12:end) = -1;
% ub(end-12:end) = 1;

lb = [];
ub = [];
nonlcon = @(C) nonl_constraint(C,r_du,r_dv,L,N_coeff,constraint_eval_loc,sensor_dB);
[C,fval,~,output] = fmincon(@(C) objective_fun(C,r_du,r_dv,r_ddu,r_dudv,r_ddv,sensor_dB',dB_per_cm_fit,wg_loc,L,N_wg),...
                                   C0,[],[],[],[],lb,ub,nonlcon,options);

% Use var C that minimize fval
if ~isempty(output.bestfeasible) && output.bestfeasible.fval < fval
    C = output.bestfeasible.x;
    fval = output.bestfeasible.fval;
end

% Keep same C0 if it does better than new C
f_ini = objective_fun(C0,r_du,r_dv,r_ddu,r_dudv,r_ddv,sensor_dB',dB_per_cm_fit,wg_loc,L,N_wg);
if fval > f_ini*1.05
    C = C0;
end

C_mat = reshape(C,3,[]);

if plot_result
    P_uv = @(A,u,v) A*[u; v; u.^2; u.*v; v.^2; u.^3; u.^2.*v; u.*v.^2; v.^3];
    evalc('colors = colormap(slanCM(''twilight''));');
    color_wg = @(dB) colors( min([max([int32(1+255*(dB+4.5)./9),1]) 255]), :);
    plot_curve_surf(C_mat,P_uv,N_wg,wg_loc,sensor_dB,view_vec,color_wg)
    cb = colorbar;caxis([-4.5 4.5]);yl = ylabel(cb,'Power (dB)','FontSize',14,'Rotation',270);yl.Position(1) = yl.Position(1)+1;
    cb.Position(1) = cb.Position(1) +0.03;
end
end

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Functions %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [f,gradf] = objective_fun(C,r_du,r_dv,r_ddu,r_dudv,r_ddv,sensor_dB,dB_per_cm_fit,wg_loc,L,N_wg)
C_mat = reshape(C,3,[]);

N_eval = 7;
ds = L/N_eval;
sensor_dB = sensor_dB(:);
f = 0;
pos = @(endpt,t) [t*(endpt(3)-endpt(1))+endpt(1) t*(endpt(4)-endpt(2))+endpt(2)];
tangent_vec = @(endpt) normalize([endpt(3)-endpt(1) endpt(4)-endpt(2)],'norm');
r_ddt = @(up,vp,ruu,ruv,rvv) up^2*ruu+2*up*vp*ruv+vp^2*rvv;
kappa_model = @(up,vp,ru,rv,ruu,ruv,rvv) norm(cross(up*ru+vp*rv,r_ddt(up,vp,ruu,ruv,rvv))) ./ (norm(up*ru+vp*rv).^3);
wg_orient = [-ones(1,N_wg) ones(1,N_wg)];

for k = 1:N_wg*2    % f eval for each waveguide
    dB_sum = 0;
    for i = linspace(0.1,0.9,N_eval)
%             u = i/N_eval-0.5/N_eval-0.5;
%             v = wg_loc(k);
        pt = pos(wg_loc(:,k),i); dir = tangent_vec(wg_loc(:,k));
        ru = r_du(C_mat,pt(1),pt(2));rv = r_dv(C_mat,pt(1),pt(2));
        ruu = r_ddu(C_mat,pt(1),pt(2));ruv = r_dudv(C_mat,pt(1),pt(2));rvv = r_ddv(C_mat,pt(1),pt(2));
        curve_direction = wg_orient(k)*sign(dot( cross(ru,rv) , r_ddt(dir(1),dir(2),ruu,ruv,rvv)) );
        optimized_kappa = curve_direction*kappa_model( dir(1),dir(2),ru,rv,ruu,ruv,rvv );

        dB_sum = dB_sum + dB_per_cm_fit(optimized_kappa*100,k)*ds;
    end
%     dB_sum
    f = f + 5*abs(sensor_dB(k) - dB_sum);
end
gradf = [];
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [c,ceq,gradc,gradceq] = nonl_constraint(C,r_du,r_dv,L,N_coeff,eval_loc,sensor_dB)
N_loc = size(eval_loc,2);
c = zeros(1,N_loc^2+4);
ceq = zeros(1,2*N_loc*N_loc);
C_mat = reshape(C,3,[]);

% no stretch constraint
for k = 1:2*N_loc
    for i = 1:N_loc
        if k <= N_loc % u1~u7
            ceq((k-1)*N_loc+i) = norm( r_du(C_mat, eval_loc(i),eval_loc(k)) ) - L;
        else    % v1~v7
            ceq((k-1)*N_loc+i) = norm( r_dv(C_mat, eval_loc(k-N_loc),eval_loc(i)) ) - L;
        end
    end
end

% ceq = ceq*10;

% no shear constraint
for i = 1:N_loc
    for j = 1:N_loc
        c((i-1)*N_loc + j) = 50* abs(dot( normalize(r_du(C_mat, eval_loc(i),eval_loc(j)),'norm'), normalize(r_dv(C_mat, eval_loc(i),eval_loc(j)),'norm') ))-1;
    end
end

% % no shear constraint
% for i = 1:N_loc
%     for j = 1:N_loc
%         ceq(2*N_wg*N_loc + (i-1)*N_loc + j) = 10*dot(r_du(C_mat, eval_loc(i),eval_loc(j)), r_dv(C_mat, eval_loc(i),eval_loc(j)));
%     end
% end

% % no saddle constraint
% c(N_loc^2+1) = -C_mat(3,3)*C_mat(3,5);
% c(N_loc^2+2) =  C_mat(3,4)^2 - 4*C_mat(3,3)*C_mat(3,5) - 1;
% c(N_loc^2+3) = -C_mat(3,4)^2 + 4*C_mat(3,3)*C_mat(3,5) - 1;
% if sum(sensor_dB(1,[2 4 6]))/3 > sum(sensor_dB(1,[1 3 5 7]))/4+0.1 && sum(sensor_dB(2,[2 4 6]))/3 > sum(sensor_dB(2,[1 3 5 7]))/4+0.1
%     c(N_loc^2+4) =  C_mat(3,3)*C_mat(3,4);
% elseif sum(sensor_dB(1,[2 4 6]))/3 < sum(sensor_dB(1,[1 3 5 7]))/4-0.1 && sum(sensor_dB(2,[2 4 6]))/3 < sum(sensor_dB(2,[1 3 5 7]))/4-0.1
%     c(N_loc^2+4) = -C_mat(3,3)*C_mat(3,4); 
% else
%     c(N_loc^2+1:N_loc^2+4) = 0;
% %     c(N_loc^2+4) = 0;
% end

% c(N_loc^2+1) = -C_mat(3,3)*C_mat(3,5);
% % c(N_loc^2+2) =  C_mat(3,3)*C_mat(3,4);
% if sum(sensor_dB(1,[2 4 6]))/3 > sum(sensor_dB(1,[1 3 5 7]))/4 && sum(sensor_dB(2,[2 4 6]))/3 < sum(sensor_dB(2,[1 3 5 7]))/4
%     c(N_loc^2+2) =  C_mat(3,4)-2*sqrt(C_mat(3,3)*C_mat(3,5))-1;
%     c(N_loc^2+3) = -C_mat(3,4)+2*sqrt(C_mat(3,3)*C_mat(3,5))-1;
% elseif sum(sensor_dB(1,[2 4 6]))/3 < sum(sensor_dB(1,[1 3 5 7]))/4 && sum(sensor_dB(2,[2 4 6]))/3 > sum(sensor_dB(2,[1 3 5 7]))/4
%     c(N_loc^2+2) =  C_mat(3,4)+2*sqrt(C_mat(3,3)*C_mat(3,5))-1;
%     c(N_loc^2+3) = -C_mat(3,4)-2*sqrt(C_mat(3,3)*C_mat(3,5))-1;
% else
%     c(N_loc^2+2) = 0;
%     c(N_loc^2+3) = 0;
% end

% gradient of ceq and c
if nargout > 2
    gradc = 50/L^2 * ones(size(C,1),size(c,2));
    for i = 1:N_loc
        for j = 1:N_loc
            u = eval_loc(i);v = eval_loc(j);
            K1 = r_du(C_mat, u,v);
            K2 = r_dv(C_mat, u,v);
            gradc(:,(i-1)*N_loc + j) = gradc(:,(i-1)*N_loc + j) .* ...
                sign(dot( K1, K2 )) .* ...
                reshape([r_dv(K1(1),u,v)'+r_du(K2(1),u,v)';r_dv(K1(2),u,v)'+r_du(K2(2),u,v)';r_dv(K1(3),u,v)'+r_du(K2(3),u,v)'],[],1);
        end
    end
    gradc(:,N_loc^2+(1:4)) = 0;
    gradc([9 12 15],N_loc^2+(1:4)) = [-C_mat(3,5) C_mat(3,4)  4*C_mat(3,5) -4*C_mat(3,5);...
                                      0           C_mat(3,3) -2*C_mat(3,4)  2*C_mat(3,4);
                                      -C_mat(3,3) 0           4*C_mat(3,3) -4*C_mat(3,3)];

    gradceq = zeros(size(C,1),2*N_loc*N_loc);
    for k = 1:2*N_loc
        for i = 1:N_loc
            if k <= N_loc
                u = eval_loc(i); v = eval_loc(k);
                Kxyz = r_du(C_mat, u,v);
                K = sum( Kxyz.^2 )^(-1/2);
                for axes = 1:3 % x,y,z
                    gradceq((0:N_coeff-1)*3+axes,(k-1)*N_loc+i) = r_du(K*Kxyz(axes),u,v);
                end
            else
                u = eval_loc(k-N_loc); v = eval_loc(i);
                Kxyz = r_dv(C_mat, u,v);
                K = sum( Kxyz.^2 )^(-1/2);
                for axes = 1:3 % x,y,z
                    gradceq((0:N_coeff-1)*3+axes,(k-1)*N_loc+i) = r_dv(K*Kxyz(axes),u,v);
                end
            end
        end
    end
end

end