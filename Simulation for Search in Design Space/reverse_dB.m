function sensor_dB = reverse_dB(C_mat,L,N_wg,wg_loc,dB_per_cm_fit)
P_uv = @(A,u,v) A*[u; v; u.^2; u.*v; v.^2; u.^3; u.^2.*v; u.*v.^2; v.^3];

r_du = @(A,u,v) A*[ones(1,size(u,2));  zeros(1,size(u,2)); 2.*u; v; zeros(1,size(u,2)); 3.*u.^2; 2.*u.*v; v.^2;    zeros(1,size(u,2))];
r_dv = @(A,u,v) A*[zeros(1,size(u,2)); ones(1,size(u,2)); zeros(1,size(u,2)); u; 2.*v;  zeros(1,size(u,2)); u.^2; 2.*u.*v; 3.*v.^2];

r_ddu = @(A,u,v) A(:,[3 6 7])*[2*ones(1,size(u,2)); 6.*u; 2.*v];
r_dudv= @(A,u,v) A(:,[4 7 8])*[ones(1,size(u,2)); 2.*u; 2.*v];
r_ddv = @(A,u,v) A(:,[5 8 9])*[2*ones(1,size(u,2)); 2.*u; 6.*v];

N_eval = 14;
ds = L/N_eval;
sensor_dB = zeros(1,2*N_wg);

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
    sensor_dB(k) = dB_sum;
end
sensor_dB = reshape(sensor_dB',[],2)';

end