% Load the state and covariance data to use to estimate the delta_v
% maneuvers
X_OD = readmatrix('X_kalman.xlsx');
P = readmatrix('covariance_history.xlsx');

% Load the time data from the observation files
filename = 'ASTE583_Project_LTB_Measurements_0-6D_Truth.csv';
f = fopen(filename, 'rt');

% Read the first line
tline = fgetl(f);

% Loop through the file line by line until the end
N = 0;
while ischar(tline)  
    % Read the next line
    tline = fgetl(f);
    N = N + 1;

    if(tline == -1)
        continue;
    end

    data_cell_array = split(tline, [","]);
    time_value = str2double(data_cell_array(1));
    time_from_file1(N) = time_value;
end
fclose(f);

% Load the time data for the rest of the days
filename2 = 'ASTE583_Project_LTB_Measurements_6D-14D_Truth.csv';
f2 = fopen(filename2, 'rt');

% Read the first line
tline2 = fgetl(f2);

% Loop through the file line by line until the end
N2 = 0;
while ischar(tline2)  
    % Read the next line
    tline2 = fgetl(f2);
    N2 = N2 + 1;

    if(tline2 == -1)
        continue;
    end

    data_cell_array2 = split(tline2, [","]);
    time_value2 = str2double(data_cell_array2(1));
    time_from_file2(N2) = time_value2;
end
fclose(f2);

% In the raw data, there are spikes at two data points that we want to
% ignore
ignored_data_points = [61, 1627];
time_from_file1(ignored_data_points) = [];
time_from_file = [time_from_file1, time_from_file2];

% Create variables and initialize state transition matrix
X0 = X_OD(end, 1:10);
n = numel(X0);
STM_0 = eye(n);
STM_0_vec = reshape(STM_0, 1, n^2);

% Set tolerance levels
tol = 1E-10;
options = odeset('RelTol',tol, 'AbsTol',tol);

% Mice path % TO BE CHANGED BY USER!
addpath('C:\Users\Brian\Desktop\College Stuff\USC\Year 2\Fall\ASTE 583 - Space Navigation Principles and Practice\mice\src\mice');
addpath('C:\Users\Brian\Desktop\College Stuff\USC\Year 2\Fall\ASTE 583 - Space Navigation Principles and Practice\mice\lib');
addpath('C:\Users\Brian\Desktop\College Stuff\USC\Year 2\Fall\ASTE 583 - Space Navigation Principles and Practice\mice\include');

% Load DE440 ephemeris % PATH TO BE CHANGED BY USER!
cspice_furnsh('C:\Users\Brian\Desktop\College Stuff\USC\Year 2\Fall\ASTE 583 - Space Navigation Principles and Practice\mice\de440s.bsp');

% Load LEAPSECOND KERNEL % PATH TO BE CHANGED BY USER!
cspice_furnsh('C:\Users\Brian\Desktop\College Stuff\USC\Year 2\Fall\ASTE 583 - Space Navigation Principles and Practice\mice\naif0012.tls.pc');

date_et = '2025 DEC 01 00:00:00 UTC' ; % start date
ET0 = cspice_str2et(date_et); % epoch time (seconds past J2000)

deg2rad = pi/180;
rad2deg = 1/deg2rad;

% Initialize additional spacecraft variables
mu_SUN = 1.3271E11;
mu_EARTH = 3.986E5;
pSRP = 4.54E-6*1E-3;
A_SC = 1.5;
m = 200;
rho_r = 0.7;
ep = (23 + 26/60 + 21.448/3600)*deg2rad;
R_EARTH = 6378;
w_EARTH = 7.292116E-5;
t0 = 0;

% Ground station location: we'll use ground station 3
lam = 40.241355*deg2rad;
phi_E = -4.2480085*deg2rad;
phiG0 = (10/60+43/3600)*deg2rad;

% Initialize the canned LTM
delta_V_LTM = [0.6931075; -0.8462091; 0.0956979]';

% Define a step size for propagation, and integrate reference trajectory
% for the first 15 days: This is before the LTM
% Propagating from FINAL OD solution (Dec 15 - Dec 16)
t_span_BEFORE_LTM = 14*86400:100:15*86400;
[~, X_BEFORE_LTM] = ode45(@EoM, t_span_BEFORE_LTM, [X0, STM_0_vec], options, ET0, mu_SUN, mu_EARTH, pSRP, A_SC, rho_r, m, n);

% Propagating using FINAL OD solution for both with/out LTM (Dec 16 - Dec 17)
t_span_POST_LTM = 15*86400:1000:16*86400;
[~, X_WITHOUT_LTM] = ode45(@EoM, t_span_POST_LTM, [X_BEFORE_LTM(end, 1:10), STM_0_vec], options, ET0, mu_SUN, mu_EARTH, pSRP, A_SC, rho_r, m, n);
[~, X_WITH_LTM] = ode45(@EoM, t_span_POST_LTM, [X_BEFORE_LTM(end, 1:10) + [0, 0, 0, delta_V_LTM, 0, 0, 0, 0], STM_0_vec], options, ET0, mu_SUN, mu_EARTH, pSRP, A_SC, rho_r, m, n);

% Obtain range and range rate measurements
for i = 1:length([t_span_BEFORE_LTM, t_span_POST_LTM])-1
    if(i < length(t_span_BEFORE_LTM))
        t(i) = t_span_BEFORE_LTM(i);
    else
        t(i) = t_span_POST_LTM(i - length(t_span_BEFORE_LTM) + 1);
    end
    phi_val = phiG0 + w_EARTH*(t(i) - t0) + phi_E;
    
    r_S_EARTH_EME(i, 1:3) = R_EARTH*[cos(phi_val)*cos(lam); sin(phi_val)*cos(lam); sin(lam)];
    v_S_EARTH_EME(i, 1:3) = R_EARTH*w_EARTH*[-sin(phi_val)*cos(lam); cos(phi_val)*cos(lam); 0];
    R_EME_to_EMO = [1, 0, 0; 0, cos(ep), sin(ep); 0, -sin(ep), cos(ep)];
    R = [R_EME_to_EMO, zeros(size(R_EME_to_EMO)); zeros(size(R_EME_to_EMO)), R_EME_to_EMO];
    X_S_EARTH_EMO(i, 1:6) = R*[r_S_EARTH_EME(i, 1:3)'; v_S_EARTH_EME(i, 1:3)'];

    r_S_EARTH_EMO(i, 1:3) = X_S_EARTH_EMO(i, 1:3);
    v_S_EARTH_EMO(i, 1:3) = X_S_EARTH_EMO(i, 4:6);

    xEarth_ephemeris = cspice_spkezr('Earth',ET0+t(i), 'ECLIPJ2000', 'NONE','Sun');
    r_13_vec(i, 1:3) = xEarth_ephemeris(1:3);
    v_13_vec(i, 1:3) = xEarth_ephemeris(4:6);

    r_S_SUN_EMO(i, 1:3) = r_S_EARTH_EMO(i, 1:3) + r_13_vec(i, 1:3);
    v_S_SUN_EMO(i, 1:3) = v_S_EARTH_EMO(i, 1:3) + v_13_vec(i, 1:3);

    if(i < length(t_span_BEFORE_LTM))
        rho_BEFORE_LTM(i) = norm([(X_BEFORE_LTM(i, 1) - r_S_SUN_EMO(i, 1)), (X_BEFORE_LTM(i, 2) - r_S_SUN_EMO(i, 2)), (X_BEFORE_LTM(i, 3) - r_S_SUN_EMO(i, 3))]);
        rho_dot_BEFORE_LTM(i) = dot([X_BEFORE_LTM(i, 1) - r_S_SUN_EMO(i, 1); X_BEFORE_LTM(i, 2) - r_S_SUN_EMO(i, 2); X_BEFORE_LTM(i, 3) - r_S_SUN_EMO(i, 3)], [X_BEFORE_LTM(i, 4) - v_S_SUN_EMO(i, 1); X_BEFORE_LTM(i, 5) - v_S_SUN_EMO(i, 2); X_BEFORE_LTM(i, 6) - v_S_SUN_EMO(i, 3)])/rho_BEFORE_LTM(i);
    else
        j = i - length(t_span_BEFORE_LTM) + 1;
        rho_WITHOUT_LTM(j) = norm([(X_WITHOUT_LTM(j, 1) - r_S_SUN_EMO(j, 1)), (X_WITHOUT_LTM(j, 2) - r_S_SUN_EMO(j, 2)), (X_WITHOUT_LTM(j, 3) - r_S_SUN_EMO(j, 3))]);
        rho_dot_WITHOUT_LTM(j) = dot([X_WITHOUT_LTM(j, 1) - r_S_SUN_EMO(j, 1); X_WITHOUT_LTM(j, 2) - r_S_SUN_EMO(j, 2); X_WITHOUT_LTM(j, 3) - r_S_SUN_EMO(j, 3)], [X_WITHOUT_LTM(j, 4) - v_S_SUN_EMO(j, 1); X_WITHOUT_LTM(j, 5) - v_S_SUN_EMO(j, 2); X_WITHOUT_LTM(j, 6) - v_S_SUN_EMO(j, 3)])/rho_WITHOUT_LTM(j);
        rho_WITH_LTM(j) = norm([(X_WITH_LTM(j, 1) - r_S_SUN_EMO(j, 1)), (X_WITH_LTM(j, 2) - r_S_SUN_EMO(j, 2)), (X_WITH_LTM(j, 3) - r_S_SUN_EMO(j, 3))]);
        rho_dot_WITH_LTM(j) = dot([X_WITH_LTM(j, 1) - r_S_SUN_EMO(j, 1); X_WITH_LTM(j, 2) - r_S_SUN_EMO(j, 2); X_WITH_LTM(j, 3) - r_S_SUN_EMO(j, 3)], [X_WITH_LTM(j, 4) - v_S_SUN_EMO(j, 1); X_WITH_LTM(j, 5) - v_S_SUN_EMO(j, 2); X_WITH_LTM(j, 6) - v_S_SUN_EMO(j, 3)])/rho_WITH_LTM(j);
    end
end

% Comute the difference between the range rate with and without the
% maneuver to evaluate the maneuver performance in real time
rho_dot_DIFF = [rho_dot_BEFORE_LTM, rho_dot_WITH_LTM] - [rho_dot_BEFORE_LTM, rho_dot_WITHOUT_LTM];
figure(1);
plot([t_span_BEFORE_LTM(2:end), t_span_POST_LTM]./86400, rho_dot_DIFF, 'bo');
title("Range Rate Residual")
xlabel("Time (days)")
ylabel("Range Rate Residual (km/s) vs time");
grid on;

% Reference Trajectory: so that we can compute the delta_v for the lunar
% correction maneuver (LCM)
R0_ref = [1.067623147085261; 1.148757045773147; -0.000321627221208]*1e8;
V0_ref = [-22.148505873534173; 18.814312217049999; -0.098774507382220];
kSRP_ref = 1;
lam_4_ref = -80*deg2rad;
phi_E4_ref = 0*deg2rad;
bias_ref = 0;
X0_ref = [R0_ref; V0_ref; kSRP_ref; lam_4_ref; phi_E4_ref; bias_ref];

% Objective: LCM to correct the estimated trajectory after 30 days (Jan 24, 2026)
% Let's do this in phases
% Integrating trajectory from Dec 1 - Dec 15
t_0_to_14 = linspace(0, 14*86400, 5000);
[~, X_0_to_14] = ode45(@EoM, t_0_to_14, [X0_ref', STM_0_vec], options, ET0, mu_SUN, mu_EARTH, pSRP, A_SC, rho_r, m, n);
% Integrating trajectory from Dec 15 - Dec 16
[~, X_14_to_15] = ode45(@EoM, t_span_BEFORE_LTM, [X_0_to_14(end, 1:n), STM_0_vec], options, ET0, mu_SUN, mu_EARTH, pSRP, A_SC, rho_r, m, n);
% Integrating trajectory from Dec 19 - Dec 25
X15_LTM = X_14_to_15(end, 1:n) + [0; 0; 0; delta_V_LTM'; 0; 0; 0; 0]';
t_span_BEFORE_LCM = linspace(15*86400, 24*86400, 4000);
[~, X_15_to_24] = ode45(@EoM, t_span_BEFORE_LCM, [X15_LTM, STM_0_vec], options, ET0, mu_SUN, mu_EARTH, pSRP, A_SC, rho_r, m, n);
% Integrating trajectory from Dec 25 - Jan 24 (30 days after LCM)
t_span_WITHOUT_LCM = 24*86400:1000:(24+30)*86400;
[~, X_24_to_54] = ode45(@EoM, t_span_WITHOUT_LCM, [X_15_to_24(end, 1:n), STM_0_vec], options, ET0, mu_SUN, mu_EARTH, pSRP, A_SC, rho_r, m, n);

% Use this state transition matrix for the LCM (going from 0-30 days)
STM = reshape(X_24_to_54(end, n+1:n+n^2), n, n);
STM_rr = STM(1:3, 1:3);
STM_rv = STM(1:3, 4:6);

% Using covariance from FINAL OD solution: Dec 15
P_day_14 = reshape(P(end, :), n, n);

% Lets determine if we meet the along-track error requirement
STM_14_to_15 = reshape(X_14_to_15(end, n+1:n+n^2), n, n);
P_15_XYZ = STM_14_to_15'*P_day_14*STM_14_to_15;

R_unit = X_14_to_15(end, 1:3)'/norm(X_14_to_15(end, 1:3)');
N_unit = cross(X_14_to_15(end, 1:3)', X_14_to_15(end, 4:6)')/norm(cross(X_14_to_15(end, 1:3)', X_14_to_15(end, 4:6)'));
T_unit = cross(N_unit, R_unit);
Rot_matrix_RTN_to_ECI = [R_unit, T_unit, N_unit];
Rot_matrix_ECI_to_RTN = inv(Rot_matrix_RTN_to_ECI);

R_matrix_ECI_to_RTN = [Rot_matrix_ECI_to_RTN, zeros(3, 3); zeros(3, 3), Rot_matrix_ECI_to_RTN];
R_matrix_ECI_to_RTN(7, 7) = 1;
R_matrix_ECI_to_RTN(8, 8) = 1;
R_matrix_ECI_to_RTN(9, 9) = 1;
R_matrix_ECI_to_RTN(10, 10) = 1;
P_RTN = R_matrix_ECI_to_RTN*P_15_XYZ*R_matrix_ECI_to_RTN';

sigma_R_15 = sqrt(P_RTN(1, 1));
sigma_T_15 = sqrt(P_RTN(2, 2));
sigma_N_15 = sqrt(P_RTN(3, 3));

% Now we're at the LTM: inplement uncertainty for sigma_delta_V_LTM
sigma_dV = 5E-3;
P_dV = zeros(n, n);
P_dV(4, 4) = sigma_dV^2;
P_dV(5, 5) = sigma_dV^2;
P_dV(6, 6) = sigma_dV^2;

% Create 1000 samples of the aposteriori covariance matrix at Dec 15 from
% the OD solution
N_samples = 1000;
X_nom = X_OD(end, 1:n)';
P_14 = reshape(P(end, :), n, n);
[U_14, D_14] = eig(P_14);

% Generate the samples at Dec 15 before the LTM happens and then propagate
X_rand = randn(n, N_samples);
X_samples = X_nom + U_14*D_14^(1/2)*X_rand; 

% Once we go from Dec 15 to Dec 16, we'll update the samples by
% implementing the sigma_dV uncertainty as noise
[U_dV, D_dV] = eig(P_dV);
dV_noise = U_dV*D_dV^(1/2)*randn(n, N_samples);

for i = 1:N_samples
    % Propagate the samples from Dec 15 to Dec 16
    [~, X_BEFORE_LTM_SAMPLE] = ode45(@EoM, t_span_BEFORE_LTM, [X_samples(:, i)', STM_0_vec], options, ET0, mu_SUN, mu_EARTH, pSRP, A_SC, rho_r, m, n);
    X_pre_LCM = X_BEFORE_LTM_SAMPLE(end, 1:10);

    % We're at Dec 16, 2025: now implement the sigma_dV uncertainty (noise)
    X_BEFORE_LCM_SAMPLE_15(i, :) = X_pre_LCM + [0, 0, 0, delta_V_LTM, 0, 0, 0, 0] + dV_noise(:, i)';

    % Propagate from Dec 16 to Dec 25 wiht the updated samples
    [~, X_BEFORE_LCM_SAMPLE] = ode45(@EoM, t_span_BEFORE_LCM, [X_BEFORE_LCM_SAMPLE_15(i, :), STM_0_vec], options, ET0, mu_SUN, mu_EARTH, pSRP, A_SC, rho_r, m, n);

    % We're at Dec 25 (the LCM): now we'll evaluate the state against the
    % state of the reference trajectory at this time
    x(i) = X_BEFORE_LCM_SAMPLE(end, 1);
    y(i) = X_BEFORE_LCM_SAMPLE(end, 2);
    z(i) = X_BEFORE_LCM_SAMPLE(end, 3);
    u(i) = X_BEFORE_LCM_SAMPLE(end, 4);
    v(i) = X_BEFORE_LCM_SAMPLE(end, 5);
    w(i) = X_BEFORE_LCM_SAMPLE(end, 6);
    r_pos(i) = norm([x(i), y(i), z(i)]);
    v_pos(i) = norm([u(i), v(i), w(i)]);

    % Calculating trajectory correction maneuver
    % First we'll calculate the difference between the estimated trajectory
    % with the reference trajectory
    R_LCM_DIFF(i, :) = [X_BEFORE_LCM_SAMPLE(end, 1) - X_15_to_24(end, 1), X_BEFORE_LCM_SAMPLE(end, 2) - X_15_to_24(end, 2), X_BEFORE_LCM_SAMPLE(end, 3) - X_15_to_24(end, 3)];
    V_LCM_DIFF(i, :) = [X_BEFORE_LCM_SAMPLE(end, 4) - X_15_to_24(end, 4), X_BEFORE_LCM_SAMPLE(end, 5) - X_15_to_24(end, 5), X_BEFORE_LCM_SAMPLE(end, 6) - X_15_to_24(end, 6)];
    
    % Use the state transtion matrix of the reference trajectory to compute
    % the required delta_V of the LCM maneuver. Also obtain its magnitude
    % and express in m/s
    delta_V_LCM(i, 1:3) = -inv(STM_rv)*STM_rr*R_LCM_DIFF(i, :)' - V_LCM_DIFF(i, :)';
    delta_V_LCM_MAG(i) = norm(delta_V_LCM(i, 1:3))*1000;

    % Let's find the total delta_V_LTM that includes its uncertainty and do
    % a monte-carlo analysis like with delta_V_LCM
    delta_V_LTM_after(i, 1:3) = delta_V_LTM + dV_noise(4:6, i)';
    delta_V_LTM_MAG(i) = norm(delta_V_LTM_after(i, 1:3));

    % Add the delta_V of the LTM and LCM to get the total delta_V for the
    % mission: this will give us a statistical analysis on the delta_V
    delta_V_total(i) = delta_V_LCM_MAG(i) + delta_V_LTM_MAG(i)*1000;
end

% Histogram plot of the spacecraft position vector difference compared to
% the position vector of the reference trajectory at Dec 25
figure(2);
histogram(r_pos - norm([X_15_to_24(end, 1), X_15_to_24(end, 2), X_15_to_24(end, 3)]));
title("Position Error Distribution")
xlabel("Position Error (km)")
ylabel("Number of Samples")
grid on;

% Histogram plot of the spacecraft velocity vector difference compared to
% the velocity vector of the reference trajectory at Dec 25
figure(3);
histogram(v_pos - norm([X_15_to_24(end, 4), X_15_to_24(end, 5), X_15_to_24(end, 6)]));
title("Velocity Error Distribution")
xlabel("Velocity Error (km/s)")
ylabel("Number of Samples")
grid on;

% Histogram of the LCM delta_V
figure(4);
histogram(delta_V_LCM_MAG);
title("LCM Distribution")
xlabel("LCM Delta-V (m/s)")
ylabel("Number of Samples")
grid on;

% Histogram of the total delta_V: LTM + LCM
figure(5);
histogram(delta_V_total);
title("Total Delta-V Distribution")
xlabel("Total Delta-V (m/s)")
ylabel("Number of Samples")
grid on;

% Compute the 99th percentile of the total delta_V to see if we're meeting
% the requirement of staying below 1139 m/s for the mission
dv99 = prctile(delta_V_total, 99);
fprintf("ΔV99: %g m/s\n", dv99)

% Function for the equations of motion of the spacecraft orbiting the sun,
% with Earth as its third body + solar radiation pressure
function dXdt = EoM(t, X, ET0, mu_SUN, mu_EARTH, pSRP, A_SC, rho_r, m, n)
    % Earth emphemeris
    xEarth_ephemeris = cspice_spkezr('Earth',ET0+t, 'ECLIPJ2000', 'NONE','Sun');
    r_13_vec = xEarth_ephemeris(1:3);
    r_31_vec = -r_13_vec;
    r_31_norm = norm(r_31_vec);

    r_12_vec = X(1:3);
    r_12_norm = norm(r_12_vec);

    r_32_vec =  r_12_vec - r_13_vec;
    r_32_norm = norm(r_32_vec);
    kSRP = X(7);

    dXdt(1:3) = X(4:6);
    dXdt(4:6) = -mu_SUN*r_12_vec/r_12_norm^3 - mu_EARTH*(r_32_vec/r_32_norm^3 - r_31_vec/r_31_norm^3) + kSRP*pSRP*A_SC*(1+rho_r)/m*r_12_vec/r_12_norm;
    dXdt(7) = 0;
    dXdt(8) = 0;
    dXdt(9) = 0;
    dXdt(10) = 0;

    x12 = r_12_vec(1);
    y12 = r_12_vec(2);
    z12 = r_12_vec(3);
    x31 = r_31_vec(1);
    y31 = r_31_vec(2);
    z31 = r_31_vec(3);

    A_jacobian = [                                                                                                                                                                                                                                                                                                                                                                                             0,                                                                                                                                                                                                                                                                                                                                                                                              0,                                                                                                                                                                                                                                                                                                                                                                                              0, 1, 0, 0,                                                             0, 0, 0, 0;...
                                                                                                                                                                                                                                                                                                                                                                                             0,                                                                                                                                                                                                                                                                                                                                                                                              0,                                                                                                                                                                                                                                                                                                                                                                                              0, 0, 1, 0,                                                             0, 0, 0, 0;...
                                                                                                                                                                                                                                                                                                                                                                                             0,                                                                                                                                                                                                                                                                                                                                                                                              0,                                                                                                                                                                                                                                                                                                                                                                                              0, 0, 0, 1,                                                             0, 0, 0, 0;...
(3*mu_SUN*x12^2)/(x12^2 + y12^2 + z12^2)^(5/2) - mu_EARTH*(1/((x12 + x31)^2 + (y12 + y31)^2 + (z12 + z31)^2)^(3/2) - (3*(2*x12 + 2*x31)*(x12 + x31))/(2*((x12 + x31)^2 + (y12 + y31)^2 + (z12 + z31)^2)^(5/2))) - mu_SUN/(x12^2 + y12^2 + z12^2)^(3/2) + (A_SC*kSRP*pSRP*(rho_r + 1))/(m*(x12^2 + y12^2 + z12^2)^(1/2)) - (A_SC*kSRP*pSRP*x12^2*(rho_r + 1))/(m*(x12^2 + y12^2 + z12^2)^(3/2)),                                                                                                                                                                 (3*mu_SUN*x12*y12)/(x12^2 + y12^2 + z12^2)^(5/2) + (3*mu_EARTH*(2*y12 + 2*y31)*(x12 + x31))/(2*((x12 + x31)^2 + (y12 + y31)^2 + (z12 + z31)^2)^(5/2)) - (A_SC*kSRP*pSRP*x12*y12*(rho_r + 1))/(m*(x12^2 + y12^2 + z12^2)^(3/2)),                                                                                                                                                                 (3*mu_EARTH*(2*z12 + 2*z31)*(x12 + x31))/(2*((x12 + x31)^2 + (y12 + y31)^2 + (z12 + z31)^2)^(5/2)) + (3*mu_SUN*x12*z12)/(x12^2 + y12^2 + z12^2)^(5/2) - (A_SC*kSRP*pSRP*x12*z12*(rho_r + 1))/(m*(x12^2 + y12^2 + z12^2)^(3/2)), 0, 0, 0, (A_SC*pSRP*x12*(rho_r + 1))/(m*(x12^2 + y12^2 + z12^2)^(1/2)), 0, 0, 0;...
                                                                                                                                                                (3*mu_SUN*x12*y12)/(x12^2 + y12^2 + z12^2)^(5/2) + (3*mu_EARTH*(2*x12 + 2*x31)*(y12 + y31))/(2*((x12 + x31)^2 + (y12 + y31)^2 + (z12 + z31)^2)^(5/2)) - (A_SC*kSRP*pSRP*x12*y12*(rho_r + 1))/(m*(x12^2 + y12^2 + z12^2)^(3/2)), (3*mu_SUN*y12^2)/(x12^2 + y12^2 + z12^2)^(5/2) - mu_EARTH*(1/((x12 + x31)^2 + (y12 + y31)^2 + (z12 + z31)^2)^(3/2) - (3*(2*y12 + 2*y31)*(y12 + y31))/(2*((x12 + x31)^2 + (y12 + y31)^2 + (z12 + z31)^2)^(5/2))) - mu_SUN/(x12^2 + y12^2 + z12^2)^(3/2) + (A_SC*kSRP*pSRP*(rho_r + 1))/(m*(x12^2 + y12^2 + z12^2)^(1/2)) - (A_SC*kSRP*pSRP*y12^2*(rho_r + 1))/(m*(x12^2 + y12^2 + z12^2)^(3/2)),                                                                                                                                                                 (3*mu_EARTH*(2*z12 + 2*z31)*(y12 + y31))/(2*((x12 + x31)^2 + (y12 + y31)^2 + (z12 + z31)^2)^(5/2)) + (3*mu_SUN*y12*z12)/(x12^2 + y12^2 + z12^2)^(5/2) - (A_SC*kSRP*pSRP*y12*z12*(rho_r + 1))/(m*(x12^2 + y12^2 + z12^2)^(3/2)), 0, 0, 0, (A_SC*pSRP*y12*(rho_r + 1))/(m*(x12^2 + y12^2 + z12^2)^(1/2)), 0, 0, 0;...
                                                                                                                                                                (3*mu_EARTH*(2*x12 + 2*x31)*(z12 + z31))/(2*((x12 + x31)^2 + (y12 + y31)^2 + (z12 + z31)^2)^(5/2)) + (3*mu_SUN*x12*z12)/(x12^2 + y12^2 + z12^2)^(5/2) - (A_SC*kSRP*pSRP*x12*z12*(rho_r + 1))/(m*(x12^2 + y12^2 + z12^2)^(3/2)),                                                                                                                                                                 (3*mu_EARTH*(2*y12 + 2*y31)*(z12 + z31))/(2*((x12 + x31)^2 + (y12 + y31)^2 + (z12 + z31)^2)^(5/2)) + (3*mu_SUN*y12*z12)/(x12^2 + y12^2 + z12^2)^(5/2) - (A_SC*kSRP*pSRP*y12*z12*(rho_r + 1))/(m*(x12^2 + y12^2 + z12^2)^(3/2)), (3*mu_SUN*z12^2)/(x12^2 + y12^2 + z12^2)^(5/2) - mu_EARTH*(1/((x12 + x31)^2 + (y12 + y31)^2 + (z12 + z31)^2)^(3/2) - (3*(2*z12 + 2*z31)*(z12 + z31))/(2*((x12 + x31)^2 + (y12 + y31)^2 + (z12 + z31)^2)^(5/2))) - mu_SUN/(x12^2 + y12^2 + z12^2)^(3/2) + (A_SC*kSRP*pSRP*(rho_r + 1))/(m*(x12^2 + y12^2 + z12^2)^(1/2)) - (A_SC*kSRP*pSRP*z12^2*(rho_r + 1))/(m*(x12^2 + y12^2 + z12^2)^(3/2)), 0, 0, 0, (A_SC*pSRP*z12*(rho_r + 1))/(m*(x12^2 + y12^2 + z12^2)^(1/2)), 0, 0, 0;...
                                                                                                                                                                                                                                                                                                                                                                                             0,                                                                                                                                                                                                                                                                                                                                                                                              0,                                                                                                                                                                                                                                                                                                                                                                                              0, 0, 0, 0,                                                             0, 0, 0, 0;...
                                                                                                                                                                                                                                                                                                                                                                                             0,                                                                                                                                                                                                                                                                                                                                                                                              0,                                                                                                                                                                                                                                                                                                                                                                                              0, 0, 0, 0,                                                             0, 0, 0, 0;...
                                                                                                                                                                                                                                                                                                                                                                                             0,                                                                                                                                                                                                                                                                                                                                                                                              0,                                                                                                                                                                                                                                                                                                                                                                                              0, 0, 0, 0,                                                             0, 0, 0, 0;...
                                                                                                                                                                                                                                                                                                                                                                                             0,                                                                                                                                                                                                                                                                                                                                                                                              0,                                                                                                                                                                                                                                                                                                                                                                                              0, 0, 0, 0,                                                             0, 0, 0, 0];
 
    dXdt(n+1:n+n^2) = A_jacobian*reshape(X(n+1:n+n^2), n, n);
    dXdt = dXdt';
end