clear; clc;
syms lam phi R_EARTH w_EARTH ep x12 y12 z12 u12 v12 w12 kSRP lam_4 phi_4 phiE4 rho_dot_bias xs ys zs us vs ws RHO phiG0 t t0 x13 y13 z13 u13 v13 w13

% Radial vector of the ground station with repect to Earth in EME2000 frame
r_s_EARTH_EME2000 = R_EARTH*[cos(phi)*cos(lam); sin(phi)*cos(lam); sin(lam)];
v_s_EARTH_EME2000 = w_EARTH*R_EARTH*[-sin(phi)*cos(lam); cos(phi)*cos(lam); 0];

% Transform from EME2000 to EMO2000
R_EME_to_EMO = [1, 0, 0; 0, cos(ep), sin(ep); 0, -sin(ep), cos(ep)];
R = [R_EME_to_EMO, zeros(size(R_EME_to_EMO)); zeros(size(R_EME_to_EMO)), R_EME_to_EMO];
X_S_EARTH_EMO = R*[r_s_EARTH_EME2000; v_s_EARTH_EME2000];

r_S_EARTH_EMO = X_S_EARTH_EMO(1:3);
v_S_EARTH_EMO = X_S_EARTH_EMO(4:6);

% Determine the coordinate of Earth with respect to the sun
r_13_vec = [x13; y13; z13];
v_13_vec = [u13; v13; w13];

r_S_SUN_EMO = r_S_EARTH_EMO + r_13_vec;
v_S_SUN_EMO = v_S_EARTH_EMO + v_13_vec;

% This is the state vector of the spacecraft we're analyzing with respect
% to the sun
r_12 = [x12; y12; z12];
v_12 = [u12; v12; w12];

rho = r_12 - r_S_SUN_EMO;
rho_mag = sqrt(rho(1)^2 + rho(2)^2 + rho(3)^2);

rho_dot = 1/rho_mag*((x12 - r_S_SUN_EMO(1))*(u12 - v_S_SUN_EMO(1)) + (y12 - r_S_SUN_EMO(2))*(v12 - v_S_SUN_EMO(2)) + (z12 - r_S_SUN_EMO(3))*(w12 - v_S_SUN_EMO(3))) + rho_dot_bias;

drho_dx = diff(rho_mag, x12);
drho_dot_dx = diff(rho_dot, x12);
%drho_dot_dx = subs(drho_dot_dx, sqrt(rho(1)^2 + rho(2)^2 + rho(3)^2), RHO);

drho_dy = diff(rho_mag, y12);
drho_dot_dy = diff(rho_dot, y12);
%drho_dot_dy = subs(drho_dot_dy, sqrt(rho(1)^2 + rho(2)^2 + rho(3)^2), RHO);

drho_dz = diff(rho_mag, z12);
drho_dot_dz = diff(rho_dot, z12);
%drho_dot_dz = subs(drho_dot_dz, sqrt(rho(1)^2 + rho(2)^2 + rho(3)^2), RHO);

drho_du = diff(rho_mag, u12);
drho_dot_du = diff(rho_dot, u12);
%drho_dot_du = subs(drho_dot_du, sqrt(rho(1)^2 + rho(2)^2 + rho(3)^2), RHO);

drho_dv = diff(rho_mag, v12);
drho_dot_dv = diff(rho_dot, v12);
%drho_dot_dv = subs(drho_dot_dv, sqrt(rho(1)^2 + rho(2)^2 + rho(3)^2), RHO);

drho_dw = diff(rho_mag, w12);
drho_dot_dw = diff(rho_dot, w12);
%drho_dot_dw = subs(drho_dot_dw, sqrt(rho(1)^2 + rho(2)^2 + rho(3)^2), RHO);

rho_updated = subs(rho_mag, [lam, phi], [lam_4, phi_4]);
rho_updated = subs(rho_updated, phi_4,  phiG0 + w_EARTH*(t - t0) + phiE4);

rho_dot_updated = subs(rho_dot, [lam, phi], [lam_4, phi_4]);
rho_dot_updated = subs(rho_dot_updated, phi_4,  phiG0 + w_EARTH*(t - t0) + phiE4);

drho_dphiE4 = diff(rho_updated, phiE4);
drho_dot_dphiE4 = diff(rho_dot_updated, phiE4);

drho_dlam4 = diff(rho_updated, lam_4);
drho_dot_dlam4 = diff(rho_dot_updated, lam_4);

drho_dkSRP = diff(rho_mag, kSRP);
drho_dot_dkSRP = diff(rho_dot, kSRP);
drho_dot_drho_dot_bias = diff(rho_dot, rho_dot_bias);
drho_drho_dot_bias = diff(rho_mag, rho_dot_bias);

H_tilde = [drho_dx, drho_dy, drho_dz, drho_du, drho_dv, drho_dw, drho_dkSRP, drho_dlam4, drho_dphiE4, drho_drho_dot_bias;...
    drho_dot_dx, drho_dot_dy, drho_dot_dz, drho_dot_du, drho_dot_dv, drho_dot_dw, drho_dot_dkSRP, drho_dot_dlam4, drho_dot_dphiE4, drho_dot_drho_dot_bias];