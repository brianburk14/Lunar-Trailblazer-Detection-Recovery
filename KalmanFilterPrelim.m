clear; close all; clc;

% Load the file for the Preliminary OD
% Read a file in text mode
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
    station_ID_value = str2double(data_cell_array(2));
    range_rate_value = str2double(data_cell_array(3));

    time_from_file(N) = time_value;
    station_ID(N) = station_ID_value;
    range_rate_from_file(N) = range_rate_value;
end
fclose(f);

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

mu_SUN = 1.3271E11;
mu_EARTH = 3.986E5;

% Define pressure in kN/m^2
pSRP = 4.54E-6*1E-3;
A_SC = 1.5;
m = 200;
R_EARTH = 6378;
w_EARTH = 7.292116E-5;
rho_r = 0.1;

lam_4_ref = -80*deg2rad;
lam_DSN = [35.244352; -35.220919; 40.241355]*deg2rad;

phi_E4_ref = 0*deg2rad;
phi_E_DSN = [-116.889538; 148.981267; -4.2480085]*deg2rad;
phiG0 = (10/60+43/3600)*deg2rad;
ep = (23 + 26/60 + 21.448/3600)*deg2rad;

% A priori values for the initial state
R0_ref = [1.067623147085261; 1.148757045773147; -0.000321627221208]*1e8;
V0_ref = [-22.148505873534173; 18.814312217049999; -0.098774507382220];
kSRP_ref = 1;
bias_ref = 0;

t0 = time_from_file(1);
X0 = [R0_ref; V0_ref; kSRP_ref; lam_4_ref; phi_E4_ref; bias_ref];
n = numel(X0);
delta_X0 = zeros(size(X0));
STM_0 = eye(n);
STM_0_vec = reshape(STM_0, n^2, 1);

% Defining known apriori covariance uncertainties for the state
% The original parameters for covariance were sigma_r = 100 km and sigma_v
% = 1E-3 km/s, but when running a passthrough with these values, we got a
% delta_X0 value with a magnitude of 1000-2000 km, which would surpass
% 10-sigma
% Let's choose our new values of 3000 km and 5E-3
sigma_r = 3000;
sigma_v = 5E-3;
sigma_kSRP = 1/3;
sigma_lam4 = 0.05*pi/180;
sigma_phiE4 = 0.05*pi/180;
sigma_bias = 1E-3;
P0_bar = diag([sigma_r^2, sigma_r^2, sigma_r^2, sigma_v^2, sigma_v^2, sigma_v^2, sigma_kSRP^2, sigma_lam4^2, sigma_phiE4^2, sigma_bias^2]);

tol = 1E-10;
options = odeset('RelTol',tol, 'AbsTol',tol);

% In the raw data, there are spikes at two data points that we want to ignore
ignored_data_points = [61, 1627];
time_from_file(ignored_data_points) = [];
station_ID(ignored_data_points) = [];
range_rate_from_file(ignored_data_points) = [];

% Downsample the data to speed up the numerical simulation
step = 15;
t = time_from_file(1:step:end);
range_rate_from_file = range_rate_from_file(1:step:end);
station_ID = station_ID(1:step:end);

% Begin looping for kalman filter
kalman_filter_iterations = 0;
criteria_met = 0;
while(criteria_met == 0)
    % Initialize the delta_X and covariance values
    delta_X_hat(:, 1) = delta_X0;
    P(:, :, 1) = P0_bar;
    X_exported(1, :) = [X0; STM_0_vec]';
    STM_final(:, :, 1) = eye(n);
    for i = 1:length(t)-1
        t_span = [t(i), t(i+1)];
        % Numerically integrate equations of motion from t(i) to t(i+1): NOT the whole trajectory!
        [~, X_out] = ode45(@EoM, t_span, [X_exported(i, 1:n)'; STM_0_vec], options, ET0, mu_SUN, mu_EARTH, pSRP, A_SC, rho_r, m, n);
        % We only care about X(i+1), which is the element at the end of X_out
        X = X_out(end, :);
        STM = reshape(X(n+1:n+n^2), n, n);
        % Lets store the estimated trajectory in X_exported
        X_exported(i+1, :) = X;

        % Obtain state transition matrix at tf so that we can map the state
        % deviation to zero later on
        STM_final(:, :, i+1) = STM*STM_final(:, :, i);

        % Propagate the apriori state deviation and covariance using
        % aposteriori information
        delta_X(:, i+1) = STM*delta_X_hat(:, i);
        P_bar(:, :, i+1) = STM*P(:, :, i)*STM';

        % Obtain locations of ground station
        switch station_ID(i+1)
            case {1, 2, 3}
                lam_val = lam_DSN(station_ID(i+1));
                phi_E_val = phi_E_DSN(station_ID(i+1));
            case 4
                lam_val = X(1, 8);
                phi_E_val = X(1, 9);
        end
        % Propagate sidereal time forward
        phi_val = phiG0 + w_EARTH*(t(i+1) - t0) + phi_E_val;
    
        % Obtain ground station position and velocity vectors in the
        % EME2000 frame, formulate the coordinate transformation between
        % EME2000 and EMO2000, and get the ground station position and
        % velocity in the EMO2000 frame
        r_S_EARTH_EME(1, 1:3) = R_EARTH*[cos(phi_val)*cos(lam_val); sin(phi_val)*cos(lam_val); sin(lam_val)];
        v_S_EARTH_EME(1, 1:3) = R_EARTH*w_EARTH*[-sin(phi_val)*cos(lam_val); cos(phi_val)*cos(lam_val); 0];
        R_EME_to_EMO = [1, 0, 0; 0, cos(ep), sin(ep); 0, -sin(ep), cos(ep)];
        R = [R_EME_to_EMO, zeros(size(R_EME_to_EMO)); zeros(size(R_EME_to_EMO)), R_EME_to_EMO];
        X_S_EARTH_EMO(1, 1:6) = R*[r_S_EARTH_EME(1, 1:3)'; v_S_EARTH_EME(1, 1:3)'];
    
        r_S_EARTH_EMO(1, 1:3) = X_S_EARTH_EMO(1, 1:3);
        v_S_EARTH_EMO(1, 1:3) = X_S_EARTH_EMO(1, 4:6);
    
        % Obtain ephemeris of Earth with respect to the sun starting from
        % Dec 1, 2025
        xEarth_ephemeris = cspice_spkezr('Earth',ET0+t(i+1), 'ECLIPJ2000', 'NONE','Sun');
        r_13_vec(1, 1:3) = xEarth_ephemeris(1:3);
        v_13_vec(1, 1:3) = xEarth_ephemeris(4:6);
    
        r_S_SUN_EMO(1, 1:3) = r_S_EARTH_EMO(1, 1:3) + r_13_vec(1, 1:3);
        v_S_SUN_EMO(1, 1:3) = v_S_EARTH_EMO(1, 1:3) + v_13_vec(1, 1:3);
    
        % Compute range and rane rate values for the model
        rho = norm([(X(1, 1) - r_S_SUN_EMO(1, 1)), (X(1, 2) - r_S_SUN_EMO(1, 2)), (X(1, 3) - r_S_SUN_EMO(1, 3))]);
        rho_dot = dot([X(1, 1) - r_S_SUN_EMO(1, 1); X(1, 2) - r_S_SUN_EMO(1, 2); X(1, 3) - r_S_SUN_EMO(1, 3)], [X(1, 4) - v_S_SUN_EMO(1, 1); X(1, 5) - v_S_SUN_EMO(1, 2); X(1, 6) - v_S_SUN_EMO(1, 3)])/rho + X(1, 10);
    
        % Obtain the model and prefit errors for the observations
        G = rho_dot;
        delta_y(1, i+1, kalman_filter_iterations+1) = range_rate_from_file(i+1) - G;
    
        % When computing the measurement mapping matrix H_tilde, we take
        % the jacobian of the observations with respect to the state in
        % calcMeasurementMappingMatrix.m as symbolic variables and define
        % those variables here
        x12 = X(1, 1);
        y12 = X(1, 2);
        z12 = X(1, 3);
        u12 = X(1, 4);
        v12 = X(1, 5);
        w12 = X(1, 6);
        x13 = r_13_vec(1, 1);
        y13 = r_13_vec(1, 2);
        z13 = r_13_vec(1, 3);
        u13 = v_13_vec(1, 1);
        v13 = v_13_vec(1, 2);
        w13 = v_13_vec(1, 3);
        lam = lam_val;
        phi = phi_val;
        drho_dot_dx = (u12 - u13 + R_EARTH*w_EARTH*cos(lam)*sin(phi))/((x13 - x12 + R_EARTH*cos(lam)*cos(phi))^2 + (y13 - y12 + R_EARTH*sin(ep)*sin(lam) + R_EARTH*cos(ep)*cos(lam)*sin(phi))^2 + (z12 - z13 - R_EARTH*cos(ep)*sin(lam) + R_EARTH*cos(lam)*sin(ep)*sin(phi))^2)^(1/2) + (((v13 - v12 + R_EARTH*w_EARTH*cos(ep)*cos(lam)*cos(phi))*(y13 - y12 + R_EARTH*sin(ep)*sin(lam) + R_EARTH*cos(ep)*cos(lam)*sin(phi)) - (x13 - x12 + R_EARTH*cos(lam)*cos(phi))*(u12 - u13 + R_EARTH*w_EARTH*cos(lam)*sin(phi)) + (w12 - w13 + R_EARTH*w_EARTH*cos(lam)*cos(phi)*sin(ep))*(z12 - z13 - R_EARTH*cos(ep)*sin(lam) + R_EARTH*cos(lam)*sin(ep)*sin(phi)))*(2*x13 - 2*x12 + 2*R_EARTH*cos(lam)*cos(phi)))/(2*((x13 - x12 + R_EARTH*cos(lam)*cos(phi))^2 + (y13 - y12 + R_EARTH*sin(ep)*sin(lam) + R_EARTH*cos(ep)*cos(lam)*sin(phi))^2 + (z12 - z13 - R_EARTH*cos(ep)*sin(lam) + R_EARTH*cos(lam)*sin(ep)*sin(phi))^2)^(3/2));
        drho_dot_dy = (((v13 - v12 + R_EARTH*w_EARTH*cos(ep)*cos(lam)*cos(phi))*(y13 - y12 + R_EARTH*sin(ep)*sin(lam) + R_EARTH*cos(ep)*cos(lam)*sin(phi)) - (x13 - x12 + R_EARTH*cos(lam)*cos(phi))*(u12 - u13 + R_EARTH*w_EARTH*cos(lam)*sin(phi)) + (w12 - w13 + R_EARTH*w_EARTH*cos(lam)*cos(phi)*sin(ep))*(z12 - z13 - R_EARTH*cos(ep)*sin(lam) + R_EARTH*cos(lam)*sin(ep)*sin(phi)))*(2*y13 - 2*y12 + 2*R_EARTH*sin(ep)*sin(lam) + 2*R_EARTH*cos(ep)*cos(lam)*sin(phi)))/(2*((x13 - x12 + R_EARTH*cos(lam)*cos(phi))^2 + (y13 - y12 + R_EARTH*sin(ep)*sin(lam) + R_EARTH*cos(ep)*cos(lam)*sin(phi))^2 + (z12 - z13 - R_EARTH*cos(ep)*sin(lam) + R_EARTH*cos(lam)*sin(ep)*sin(phi))^2)^(3/2)) - (v13 - v12 + R_EARTH*w_EARTH*cos(ep)*cos(lam)*cos(phi))/((x13 - x12 + R_EARTH*cos(lam)*cos(phi))^2 + (y13 - y12 + R_EARTH*sin(ep)*sin(lam) + R_EARTH*cos(ep)*cos(lam)*sin(phi))^2 + (z12 - z13 - R_EARTH*cos(ep)*sin(lam) + R_EARTH*cos(lam)*sin(ep)*sin(phi))^2)^(1/2);
        drho_dot_dz = (w12 - w13 + R_EARTH*w_EARTH*cos(lam)*cos(phi)*sin(ep))/((x13 - x12 + R_EARTH*cos(lam)*cos(phi))^2 + (y13 - y12 + R_EARTH*sin(ep)*sin(lam) + R_EARTH*cos(ep)*cos(lam)*sin(phi))^2 + (z12 - z13 - R_EARTH*cos(ep)*sin(lam) + R_EARTH*cos(lam)*sin(ep)*sin(phi))^2)^(1/2) - (((v13 - v12 + R_EARTH*w_EARTH*cos(ep)*cos(lam)*cos(phi))*(y13 - y12 + R_EARTH*sin(ep)*sin(lam) + R_EARTH*cos(ep)*cos(lam)*sin(phi)) - (x13 - x12 + R_EARTH*cos(lam)*cos(phi))*(u12 - u13 + R_EARTH*w_EARTH*cos(lam)*sin(phi)) + (w12 - w13 + R_EARTH*w_EARTH*cos(lam)*cos(phi)*sin(ep))*(z12 - z13 - R_EARTH*cos(ep)*sin(lam) + R_EARTH*cos(lam)*sin(ep)*sin(phi)))*(2*z12 - 2*z13 - 2*R_EARTH*cos(ep)*sin(lam) + 2*R_EARTH*cos(lam)*sin(ep)*sin(phi)))/(2*((x13 - x12 + R_EARTH*cos(lam)*cos(phi))^2 + (y13 - y12 + R_EARTH*sin(ep)*sin(lam) + R_EARTH*cos(ep)*cos(lam)*sin(phi))^2 + (z12 - z13 - R_EARTH*cos(ep)*sin(lam) + R_EARTH*cos(lam)*sin(ep)*sin(phi))^2)^(3/2));
        drho_dot_du = -(x13 - x12 + R_EARTH*cos(lam)*cos(phi))/((x13 - x12 + R_EARTH*cos(lam)*cos(phi))^2 + (y13 - y12 + R_EARTH*sin(ep)*sin(lam) + R_EARTH*cos(ep)*cos(lam)*sin(phi))^2 + (z12 - z13 - R_EARTH*cos(ep)*sin(lam) + R_EARTH*cos(lam)*sin(ep)*sin(phi))^2)^(1/2);
        drho_dot_dv = -(y13 - y12 + R_EARTH*sin(ep)*sin(lam) + R_EARTH*cos(ep)*cos(lam)*sin(phi))/((x13 - x12 + R_EARTH*cos(lam)*cos(phi))^2 + (y13 - y12 + R_EARTH*sin(ep)*sin(lam) + R_EARTH*cos(ep)*cos(lam)*sin(phi))^2 + (z12 - z13 - R_EARTH*cos(ep)*sin(lam) + R_EARTH*cos(lam)*sin(ep)*sin(phi))^2)^(1/2);
        drho_dot_dw = (z12 - z13 - R_EARTH*cos(ep)*sin(lam) + R_EARTH*cos(lam)*sin(ep)*sin(phi))/((x13 - x12 + R_EARTH*cos(lam)*cos(phi))^2 + (y13 - y12 + R_EARTH*sin(ep)*sin(lam) + R_EARTH*cos(ep)*cos(lam)*sin(phi))^2 + (z12 - z13 - R_EARTH*cos(ep)*sin(lam) + R_EARTH*cos(lam)*sin(ep)*sin(phi))^2)^(1/2);
        drho_dot_dkSRP = 0;
        drho_dot_dbias = 1;
    
        % We need to consider when ground station 4 is making the
        % measurements since it's part of the state
        if(station_ID(i+1) == 4)
            lam_4 = lam_val;
            phiE4 = phi_E_val;
            drho_dot_dlam4 = ((R_EARTH*cos(lam_4)*sin(ep) - R_EARTH*sin(phiG0 + phiE4 + w_EARTH*(t(i+1) - t0))*cos(ep)*sin(lam_4))*(v13 - v12 + R_EARTH*w_EARTH*cos(phiG0 + phiE4 + w_EARTH*(t(i+1) - t0))*cos(ep)*cos(lam_4)) - (R_EARTH*cos(ep)*cos(lam_4) + R_EARTH*sin(phiG0 + phiE4 + w_EARTH*(t(i+1) - t0))*sin(ep)*sin(lam_4))*(w12 - w13 + R_EARTH*w_EARTH*cos(phiG0 + phiE4 + w_EARTH*(t(i+1) - t0))*cos(lam_4)*sin(ep)) + R_EARTH*cos(phiG0 + phiE4 + w_EARTH*(t(i+1) - t0))*sin(lam_4)*(u12 - u13 + R_EARTH*w_EARTH*sin(phiG0 + phiE4 + w_EARTH*(t(i+1) - t0))*cos(lam_4)) + R_EARTH*w_EARTH*sin(phiG0 + phiE4 + w_EARTH*(t(i+1) - t0))*sin(lam_4)*(x13 - x12 + R_EARTH*cos(phiG0 + phiE4 + w_EARTH*(t(i+1) - t0))*cos(lam_4)) - R_EARTH*w_EARTH*cos(phiG0 + phiE4 + w_EARTH*(t(i+1) - t0))*cos(ep)*sin(lam_4)*(y13 - y12 + R_EARTH*sin(ep)*sin(lam_4) + R_EARTH*sin(phiG0 + phiE4 + w_EARTH*(t(i+1) - t0))*cos(ep)*cos(lam_4)) - R_EARTH*w_EARTH*cos(phiG0 + phiE4 + w_EARTH*(t(i+1) - t0))*sin(ep)*sin(lam_4)*(z12 - z13 - R_EARTH*cos(ep)*sin(lam_4) + R_EARTH*sin(phiG0 + phiE4 + w_EARTH*(t(i+1) - t0))*cos(lam_4)*sin(ep)))/((x13 - x12 + R_EARTH*cos(phiG0 + phiE4 + w_EARTH*(t(i+1) - t0))*cos(lam_4))^2 + (y13 - y12 + R_EARTH*sin(ep)*sin(lam_4) + R_EARTH*sin(phiG0 + phiE4 + w_EARTH*(t(i+1) - t0))*cos(ep)*cos(lam_4))^2 + (z12 - z13 - R_EARTH*cos(ep)*sin(lam_4) + R_EARTH*sin(phiG0 + phiE4 + w_EARTH*(t(i+1) - t0))*cos(lam_4)*sin(ep))^2)^(1/2) + (((v13 - v12 + R_EARTH*w_EARTH*cos(phiG0 + phiE4 + w_EARTH*(t(i+1) - t0))*cos(ep)*cos(lam_4))*(y13 - y12 + R_EARTH*sin(ep)*sin(lam_4) + R_EARTH*sin(phiG0 + phiE4 + w_EARTH*(t(i+1) - t0))*cos(ep)*cos(lam_4)) - (x13 - x12 + R_EARTH*cos(phiG0 + phiE4 + w_EARTH*(t(i+1) - t0))*cos(lam_4))*(u12 - u13 + R_EARTH*w_EARTH*sin(phiG0 + phiE4 + w_EARTH*(t(i+1) - t0))*cos(lam_4)) + (w12 - w13 + R_EARTH*w_EARTH*cos(phiG0 + phiE4 + w_EARTH*(t(i+1) - t0))*cos(lam_4)*sin(ep))*(z12 - z13 - R_EARTH*cos(ep)*sin(lam_4) + R_EARTH*sin(phiG0 + phiE4 + w_EARTH*(t(i+1) - t0))*cos(lam_4)*sin(ep)))*(2*(R_EARTH*cos(ep)*cos(lam_4) + R_EARTH*sin(phiG0 + phiE4 + w_EARTH*(t(i+1) - t0))*sin(ep)*sin(lam_4))*(z12 - z13 - R_EARTH*cos(ep)*sin(lam_4) + R_EARTH*sin(phiG0 + phiE4 + w_EARTH*(t(i+1) - t0))*cos(lam_4)*sin(ep)) - 2*(R_EARTH*cos(lam_4)*sin(ep) - R_EARTH*sin(phiG0 + phiE4 + w_EARTH*(t(i+1) - t0))*cos(ep)*sin(lam_4))*(y13 - y12 + R_EARTH*sin(ep)*sin(lam_4) + R_EARTH*sin(phiG0 + phiE4 + w_EARTH*(t(i+1) - t0))*cos(ep)*cos(lam_4)) + 2*R_EARTH*cos(phiG0 + phiE4 + w_EARTH*(t(i+1) - t0))*sin(lam_4)*(x13 - x12 + R_EARTH*cos(phiG0 + phiE4 + w_EARTH*(t(i+1) - t0))*cos(lam_4))))/(2*((x13 - x12 + R_EARTH*cos(phiG0 + phiE4 + w_EARTH*(t(i+1) - t0))*cos(lam_4))^2 + (y13 - y12 + R_EARTH*sin(ep)*sin(lam_4) + R_EARTH*sin(phiG0 + phiE4 + w_EARTH*(t(i+1) - t0))*cos(ep)*cos(lam_4))^2 + (z12 - z13 - R_EARTH*cos(ep)*sin(lam_4) + R_EARTH*sin(phiG0 + phiE4 + w_EARTH*(t(i+1) - t0))*cos(lam_4)*sin(ep))^2)^(3/2));
            drho_dot_dphiE4 = (R_EARTH*sin(phiG0 + phiE4 + w_EARTH*(t(i+1) - t0))*cos(lam_4)*(u12 - u13 + R_EARTH*w_EARTH*sin(phiG0 + phiE4 + w_EARTH*(t(i+1) - t0))*cos(lam_4)) + R_EARTH*cos(phiG0 + phiE4 + w_EARTH*(t(i+1) - t0))*cos(lam_4)*sin(ep)*(w12 - w13 + R_EARTH*w_EARTH*cos(phiG0 + phiE4 + w_EARTH*(t(i+1) - t0))*cos(lam_4)*sin(ep)) - R_EARTH*w_EARTH*cos(phiG0 + phiE4 + w_EARTH*(t(i+1) - t0))*cos(lam_4)*(x13 - x12 + R_EARTH*cos(phiG0 + phiE4 + w_EARTH*(t(i+1) - t0))*cos(lam_4)) + R_EARTH*cos(phiG0 + phiE4 + w_EARTH*(t(i+1) - t0))*cos(ep)*cos(lam_4)*(v13 - v12 + R_EARTH*w_EARTH*cos(phiG0 + phiE4 + w_EARTH*(t(i+1) - t0))*cos(ep)*cos(lam_4)) - R_EARTH*w_EARTH*sin(phiG0 + phiE4 + w_EARTH*(t(i+1) - t0))*cos(ep)*cos(lam_4)*(y13 - y12 + R_EARTH*sin(ep)*sin(lam_4) + R_EARTH*sin(phiG0 + phiE4 + w_EARTH*(t(i+1) - t0))*cos(ep)*cos(lam_4)) - R_EARTH*w_EARTH*sin(phiG0 + phiE4 + w_EARTH*(t(i+1) - t0))*cos(lam_4)*sin(ep)*(z12 - z13 - R_EARTH*cos(ep)*sin(lam_4) + R_EARTH*sin(phiG0 + phiE4 + w_EARTH*(t(i+1) - t0))*cos(lam_4)*sin(ep)))/((x13 - x12 + R_EARTH*cos(phiG0 + phiE4 + w_EARTH*(t(i+1) - t0))*cos(lam_4))^2 + (y13 - y12 + R_EARTH*sin(ep)*sin(lam_4) + R_EARTH*sin(phiG0 + phiE4 + w_EARTH*(t(i+1) - t0))*cos(ep)*cos(lam_4))^2 + (z12 - z13 - R_EARTH*cos(ep)*sin(lam_4) + R_EARTH*sin(phiG0 + phiE4 + w_EARTH*(t(i+1) - t0))*cos(lam_4)*sin(ep))^2)^(1/2) - (((v13 - v12 + R_EARTH*w_EARTH*cos(phiG0 + phiE4 + w_EARTH*(t(i+1) - t0))*cos(ep)*cos(lam_4))*(y13 - y12 + R_EARTH*sin(ep)*sin(lam_4) + R_EARTH*sin(phiG0 + phiE4 + w_EARTH*(t(i+1) - t0))*cos(ep)*cos(lam_4)) - (x13 - x12 + R_EARTH*cos(phiG0 + phiE4 + w_EARTH*(t(i+1) - t0))*cos(lam_4))*(u12 - u13 + R_EARTH*w_EARTH*sin(phiG0 + phiE4 + w_EARTH*(t(i+1) - t0))*cos(lam_4)) + (w12 - w13 + R_EARTH*w_EARTH*cos(phiG0 + phiE4 + w_EARTH*(t(i+1) - t0))*cos(lam_4)*sin(ep))*(z12 - z13 - R_EARTH*cos(ep)*sin(lam_4) + R_EARTH*sin(phiG0 + phiE4 + w_EARTH*(t(i+1) - t0))*cos(lam_4)*sin(ep)))*(2*R_EARTH*cos(phiG0 + phiE4 + w_EARTH*(t(i+1) - t0))*cos(ep)*cos(lam_4)*(y13 - y12 + R_EARTH*sin(ep)*sin(lam_4) + R_EARTH*sin(phiG0 + phiE4 + w_EARTH*(t(i+1) - t0))*cos(ep)*cos(lam_4)) - 2*R_EARTH*sin(phiG0 + phiE4 + w_EARTH*(t(i+1) - t0))*cos(lam_4)*(x13 - x12 + R_EARTH*cos(phiG0 + phiE4 + w_EARTH*(t(i+1) - t0))*cos(lam_4)) + 2*R_EARTH*cos(phiG0 + phiE4 + w_EARTH*(t(i+1) - t0))*cos(lam_4)*sin(ep)*(z12 - z13 - R_EARTH*cos(ep)*sin(lam_4) + R_EARTH*sin(phiG0 + phiE4 + w_EARTH*(t(i+1) - t0))*cos(lam_4)*sin(ep))))/(2*((x13 - x12 + R_EARTH*cos(phiG0 + phiE4 + w_EARTH*(t(i+1) - t0))*cos(lam_4))^2 + (y13 - y12 + R_EARTH*sin(ep)*sin(lam_4) + R_EARTH*sin(phiG0 + phiE4 + w_EARTH*(t(i+1) - t0))*cos(ep)*cos(lam_4))^2 + (z12 - z13 - R_EARTH*cos(ep)*sin(lam_4) + R_EARTH*sin(phiG0 + phiE4 + w_EARTH*(t(i+1) - t0))*cos(lam_4)*sin(ep))^2)^(3/2));
        else
            drho_dot_dlam4 = 0;
            drho_dot_dphiE4 = 0;
        end
    
        % Obtain measurement mapping matrix
        H_tilde(1, 1:n, i+1) = [drho_dot_dx, drho_dot_dy, drho_dot_dz, drho_dot_du, drho_dot_dv, drho_dot_dw, drho_dot_dkSRP, drho_dot_dlam4, drho_dot_dphiE4, drho_dot_dbias];
        switch station_ID(i+1)
            case {1, 2, 3}
                sigma_rho = 1E-3;
                sigma_rho_dot = 0.1E-6;
            case 4
                sigma_rho = 10E-3;
                sigma_rho_dot = 1E-6;
        end

        % Initialize measurement covariance, compute the Kalman gain, and
        % obtain the aposteriori state deviation at i+1, same goes for
        % aposteriori covariance
        R = sigma_rho_dot^2;
        K(:, :, i+1) = P_bar(:, :, i+1)*transpose(H_tilde(1, 1:n, i+1))*inv(H_tilde(1, 1:n, i+1)*P_bar(:, :, i+1)*transpose(H_tilde(1, 1:n, i+1)) + R);
        delta_X_hat(:, i+1) = delta_X(:, i+1) + K(:, :, i+1)*(delta_y(1, i+1, kalman_filter_iterations+1) - H_tilde(1, 1:n, i+1)*delta_X(:, i+1));
        P(:, :, i+1) = (eye(size(K(:, :, i+1)*H_tilde(1, 1:n, i+1))) - K(:, :, i+1)*H_tilde(1, 1:n, i+1))*P_bar(:, :, i+1)*(eye(size(K(:, :, i+1)*H_tilde(1, 1:n, i+1))) - K(:, :, i+1)*H_tilde(1, 1:n, i+1))' + K(:, :, i+1)*R*K(:, :, i+1)';

        % Compute postfit error to see how it behaves with time for each
        % iteration of the EKF
        postfit_error(1, i+1, kalman_filter_iterations+1) = delta_y(1, i+1, kalman_filter_iterations+1) - H_tilde(1, 1:n, i+1)*delta_X(:, i+1);

        sigma_x(i) = sqrt(P(1, 1, i+1));
        sigma_y(i) = sqrt(P(2, 2, i+1));
        sigma_z(i) = sqrt(P(3, 3, i+1));
        sigma_kSRP(i) = sqrt(P(7, 7, i+1));

        % Transform the covariance from the EME2000 (ECI) frame to RTN frame
        R_unit = X(1, 1:3)'/norm(X(1, 1:3)');
        N_unit = cross(X(1, 1:3)', X(1, 4:6)')/norm(cross(X(1, 1:3)', X(1, 4:6)'));
        T_unit = cross(N_unit, R_unit);
        Rot_matrix_RTN_to_ECI = [R_unit, T_unit, N_unit];
        Rot_matrix_ECI_to_RTN = inv(Rot_matrix_RTN_to_ECI);
    
        R_matrix_ECI_to_RTN = [Rot_matrix_ECI_to_RTN, zeros(3, 3); zeros(3, 3), Rot_matrix_ECI_to_RTN];
        R_matrix_ECI_to_RTN(7, 7) = 1;
        R_matrix_ECI_to_RTN(8, 8) = 1;
        R_matrix_ECI_to_RTN(9, 9) = 1;
        R_matrix_ECI_to_RTN(10, 10) = 1;
        P_RTN = R_matrix_ECI_to_RTN*P(:, :, i+1)*R_matrix_ECI_to_RTN';

        % Obtain position covariance in the RTN frame
        sigma_R(i) = sqrt(P_RTN(1, 1));
        sigma_T(i) = sqrt(P_RTN(2, 2));
        sigma_N(i) = sqrt(P_RTN(3, 3));
    end

    % Map the state deviation back to zero to compute the estimated state
    % deviation at t = 0. Update X0 and delta_X0 and see if we have converged
    delta_X0_hat = inv(STM_final(:, :, end))*delta_X_hat(:, end);
    X0 = X0 + delta_X0_hat;
    delta_X0 = delta_X0 - delta_X0_hat;

    delta_X0_hat'*inv(P0_bar)*delta_X0_hat
    criterion(kalman_filter_iterations+1) = delta_X0_hat'*inv(P0_bar)*delta_X0_hat;
    if(kalman_filter_iterations > 0)
        if(criterion(kalman_filter_iterations+1) > criterion(kalman_filter_iterations))
            error("Failed to converge");
        end
    end
    if(delta_X0_hat'*inv(P0_bar)*delta_X0_hat < tol)
        criteria_met = 1;
        disp("Converged with " + kalman_filter_iterations + " iterations.");
    else
        kalman_filter_iterations = kalman_filter_iterations + 1;
    end
end

% Some values of the postfit error have data points that are spikes, so
%those points will be ignored
ignored_time_steps = [2, 51];
t(ignored_time_steps) = [];
delta_y(:, ignored_time_steps, :) = [];
postfit_error(:, ignored_time_steps, :) = [];
sigma_R(ignored_time_steps) = [];
sigma_T(ignored_time_steps) = [];
sigma_N(ignored_time_steps) = [];

figure(1);
for i = 1:kalman_filter_iterations+1
    subplot(kalman_filter_iterations+1, 2, 2*i-1);
    plot(t(2:end)/86400, delta_y(1, 2:end, i), 'bo');
    title("Iteration " + (i-1) + ": Range Rate Prefit");
    xlabel("time (days)");
    ylabel('$\delta \dot{\rho}$ (km/s)', 'Interpreter','latex');
    grid on;
    
    subplot(kalman_filter_iterations+1, 2, 2*i);
    plot(t(2:end)/86400, postfit_error(1, 2:end, i), 'bo');
    title("Iteration " + (i-1) + ": Range Rate Postfit");
    xlabel("time (days)");
    ylabel('$\delta \dot{\rho}$ (km/s)', 'Interpreter','latex')
    grid on;
end

figure(2);
plot(t(2:end)/86400, sigma_R, 'go');
title("Spacecraft Uncertainty vs. Time")
xlabel("Time (days)")
ylabel("Error (km) vs time");
grid on;
hold on;
plot(t(2:end)/86400, sigma_T, 'ro');
hold on;
plot(t(2:end)/86400, sigma_N, 'bo');
legend("σ_R", "σ_T", "σ_N");

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