% Define symbolic variables
clear; clc;
syms mu_SUN mu_EARTH x12 y12 z12 u12 v12 w12 x31 y31 z31 r_12_norm r_32_norm r_31_norm pSRP rho_r m A_SC kSRP lam4 phiE4 rho_dot_bias;

% Define the state of the spacecraft
X = [x12; y12; z12; u12; v12; w12; kSRP; lam4; phiE4; rho_dot_bias];
dXdt(1:3) = X(4:6);
dXdt(4) = -mu_SUN*x12/r_12_norm^3 - mu_EARTH*((x12 + x31)/r_32_norm^3 + x31/r_31_norm^3) + kSRP*pSRP*(1+rho_r)*A_SC*x12/(m*r_12_norm);
dXdt(5) = -mu_SUN*y12/r_12_norm^3 - mu_EARTH*((y12 + y31)/r_32_norm^3 + y31/r_31_norm^3) + kSRP*pSRP*(1+rho_r)*A_SC*y12/(m*r_12_norm);
dXdt(6) = -mu_SUN*z12/r_12_norm^3 - mu_EARTH*((z12 + z31)/r_32_norm^3 + z31/r_31_norm^3) + kSRP*pSRP*(1+rho_r)*A_SC*z12/(m*r_12_norm);
dXdt(7) = 0;
dXdt(8) = 0;
dXdt(9) = 0;
dXdt(10) = 0;

% We'll substitute r with the sqrt(x^2 + y^2 + z^2) so that we can take the
% partial derivatives correctly
dXdt = subs(dXdt, [r_12_norm, r_31_norm, r_32_norm], [(x12^2 + y12^2 + z12^2)^(1/2),...
            (x31^2 + y31^2 + z31^2)^(1/2),...
            ((x12 + x31)^2 + (y12 + y31)^2 + (z12 + z31)^2)^(1/2)]);
for i = 1:length(X)
    for j = 1:length(X)
        A_mat(i, j) = diff(dXdt(i), X(j));
    end
end

% Checking (J is the sam as A: verified with text compare)
f = transpose(dXdt);
J = jacobian(f, X);