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
rho_r = 0.7;

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

sigma_r = 100;
sigma_v = 1E-3;
sigma_kSRP = 1/3;
sigma_lam4 = 0.001*pi/180;
sigma_phiE4 = 0.001*pi/180;
sigma_bias = 1E-4;
P0_bar = diag([sigma_r^2, sigma_r^2, sigma_r^2, sigma_v^2, sigma_v^2, sigma_v^2, sigma_kSRP^2, sigma_lam4^2, sigma_phiE4^2, sigma_bias^2]);

tol = 1E-10;
options = odeset('RelTol',tol, 'AbsTol',tol);

% In the raw data, there are spikes at two data points that we want to
% ignore
ignored_data_points = [61, 1627];
time_from_file(ignored_data_points) = [];
station_ID(ignored_data_points) = [];
range_rate_from_file(ignored_data_points) = [];

t = time_from_file(1:100:end);
range_rate_from_file = range_rate_from_file(1:100:end);
station_ID = station_ID(1:100:end);

t_span = [t(1), t(end)];
[~, X_out] = ode45(@EoM, t_span, [X0; STM_0_vec], options, ET0, mu_SUN, mu_EARTH, pSRP, A_SC, rho_r, m, n);

STM_final = reshape(X_out(end, n+1:n+n^2), n, n);
STM_prelim_Kalman_iter0 = [1.000941393662466	0.006844382026995	-3.723755807708114e-06	5.183199389400507e+05	1.186367458161943e+03	-0.569121796312836	-4.921031166844370	0	0	0
0.006843123502727	1.003734858942944	-4.598443103145806e-06	1.186259030135554e+03	5.187246925401854e+05	-0.680342597042378	-6.025308844467725	0	0	0
-3.720643498335707e-06	-4.595473502623874e-06	0.995345252274643	-0.568853655451469	-0.680086751787874	5.173175858668199e+05	0.003288244430977	0	0	0
4.499811904258522e-09	2.653244224403266e-08	-1.273176630456809e-11	1.001373973786008	0.006890111513739	-2.862877789430273e-06	-1.939284937853297e-05	0	0	0
2.652051197123343e-08	1.358167198626856e-08	-1.530767684210767e-11	0.006888890158946	1.003282349236676	-3.322530419517379e-06	-2.298388331359687e-05	0	0	0
-1.270226292239517e-11	-1.527952633636509e-11	-1.791568390023424e-08	-2.859857395885766e-06	-3.319648531717404e-06	0.995365074965944	1.107740293931128e-08	0	0	0
0	0	0	0	0	0	1	0	0	0
0	0	0	0	0	0	0	1	0	0
0	0	0	0	0	0	0	0	1	0
0	0	0	0	0	0	0	0	0	1];

STM_dev = abs(STM_final - STM_prelim_Kalman_iter0)

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
    dXdt(4:6) = -mu_SUN*r_12_vec/r_12_norm^3 - mu_EARTH*(r_32_vec/r_32_norm^3 - r_31_vec/r_31_norm^3) - kSRP*pSRP*A_SC*(1+rho_r)/m*r_12_vec/r_12_norm;
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

    A_jacobian = [0, 0, 0, 1, 0, 0, 0, 0, 0, 0;...
            0, 0, 0, 0, 1, 0, 0, 0, 0, 0;...
            0, 0, 0, 0, 0, 1, 0, 0, 0, 0;...
            (3*mu_SUN*x12^2)/(x12^2 + y12^2 + z12^2)^(5/2) - mu_EARTH*(1/((x12 + x31)^2 + (y12 + y31)^2 + (z12 + z31)^2)^(3/2) - (3*(2*x12 + 2*x31)*(x12 + x31))/(2*((x12 + x31)^2 + (y12 + y31)^2 + (z12 + z31)^2)^(5/2))) - mu_SUN/(x12^2 + y12^2 + z12^2)^(3/2) - (A_SC*kSRP*pSRP*(rho_r + 1))/(m*(x12^2 + y12^2 + z12^2)^(1/2)) + (A_SC*kSRP*pSRP*x12^2*(rho_r + 1))/(m*(x12^2 + y12^2 + z12^2)^(3/2)), (3*mu_SUN*x12*y12)/(x12^2 + y12^2 + z12^2)^(5/2) + (3*mu_EARTH*(2*y12 + 2*y31)*(x12 + x31))/(2*((x12 + x31)^2 + (y12 + y31)^2 + (z12 + z31)^2)^(5/2)) + (A_SC*kSRP*pSRP*x12*y12*(rho_r + 1))/(m*(x12^2 + y12^2 + z12^2)^(3/2)), (3*mu_EARTH*(2*z12 + 2*z31)*(x12 + x31))/(2*((x12 + x31)^2 + (y12 + y31)^2 + (z12 + z31)^2)^(5/2)) + (3*mu_SUN*x12*z12)/(x12^2 + y12^2 + z12^2)^(5/2) + (A_SC*kSRP*pSRP*x12*z12*(rho_r + 1))/(m*(x12^2 + y12^2 + z12^2)^(3/2)), 0, 0, 0, -(A_SC*pSRP*x12*(rho_r + 1))/(m*(x12^2 + y12^2 + z12^2)^(1/2)), 0, 0, 0;...
            (3*mu_SUN*x12*y12)/(x12^2 + y12^2 + z12^2)^(5/2) + (3*mu_EARTH*(2*x12 + 2*x31)*(y12 + y31))/(2*((x12 + x31)^2 + (y12 + y31)^2 + (z12 + z31)^2)^(5/2)) + (A_SC*kSRP*pSRP*x12*y12*(rho_r + 1))/(m*(x12^2 + y12^2 + z12^2)^(3/2)), (3*mu_SUN*y12^2)/(x12^2 + y12^2 + z12^2)^(5/2) - mu_EARTH*(1/((x12 + x31)^2 + (y12 + y31)^2 + (z12 + z31)^2)^(3/2) - (3*(2*y12 + 2*y31)*(y12 + y31))/(2*((x12 + x31)^2 + (y12 + y31)^2 + (z12 + z31)^2)^(5/2))) - mu_SUN/(x12^2 + y12^2 + z12^2)^(3/2) - (A_SC*kSRP*pSRP*(rho_r + 1))/(m*(x12^2 + y12^2 + z12^2)^(1/2)) + (A_SC*kSRP*pSRP*y12^2*(rho_r + 1))/(m*(x12^2 + y12^2 + z12^2)^(3/2)), (3*mu_EARTH*(2*z12 + 2*z31)*(y12 + y31))/(2*((x12 + x31)^2 + (y12 + y31)^2 + (z12 + z31)^2)^(5/2)) + (3*mu_SUN*y12*z12)/(x12^2 + y12^2 + z12^2)^(5/2) + (A_SC*kSRP*pSRP*y12*z12*(rho_r + 1))/(m*(x12^2 + y12^2 + z12^2)^(3/2)), 0, 0, 0, -(A_SC*pSRP*y12*(rho_r + 1))/(m*(x12^2 + y12^2 + z12^2)^(1/2)), 0, 0, 0;...
            (3*mu_EARTH*(2*x12 + 2*x31)*(z12 + z31))/(2*((x12 + x31)^2 + (y12 + y31)^2 + (z12 + z31)^2)^(5/2)) + (3*mu_SUN*x12*z12)/(x12^2 + y12^2 + z12^2)^(5/2) + (A_SC*kSRP*pSRP*x12*z12*(rho_r + 1))/(m*(x12^2 + y12^2 + z12^2)^(3/2)), (3*mu_EARTH*(2*y12 + 2*y31)*(z12 + z31))/(2*((x12 + x31)^2 + (y12 + y31)^2 + (z12 + z31)^2)^(5/2)) + (3*mu_SUN*y12*z12)/(x12^2 + y12^2 + z12^2)^(5/2) + (A_SC*kSRP*pSRP*y12*z12*(rho_r + 1))/(m*(x12^2 + y12^2 + z12^2)^(3/2)), (3*mu_SUN*z12^2)/(x12^2 + y12^2 + z12^2)^(5/2) - mu_EARTH*(1/((x12 + x31)^2 + (y12 + y31)^2 + (z12 + z31)^2)^(3/2) - (3*(2*z12 + 2*z31)*(z12 + z31))/(2*((x12 + x31)^2 + (y12 + y31)^2 + (z12 + z31)^2)^(5/2))) - mu_SUN/(x12^2 + y12^2 + z12^2)^(3/2) - (A_SC*kSRP*pSRP*(rho_r + 1))/(m*(x12^2 + y12^2 + z12^2)^(1/2)) + (A_SC*kSRP*pSRP*z12^2*(rho_r + 1))/(m*(x12^2 + y12^2 + z12^2)^(3/2)), 0, 0, 0, -(A_SC*pSRP*z12*(rho_r + 1))/(m*(x12^2 + y12^2 + z12^2)^(1/2)), 0, 0, 0;...
            0, 0, 0, 0, 0, 0, 0, 0, 0, 0;...
            0, 0, 0, 0, 0, 0, 0, 0, 0, 0;...
            0, 0, 0, 0, 0, 0, 0, 0, 0, 0;...
            0, 0, 0, 0, 0, 0, 0, 0, 0, 0];
    dXdt(n+1:n+n^2) = A_jacobian*reshape(X(n+1:n+n^2), n, n);
    dXdt = dXdt';
end