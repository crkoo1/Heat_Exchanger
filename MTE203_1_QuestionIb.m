%-------------------------------------------------------------------------%
% Part 1: Shell and Tube Heat Exchanger Volume
% Code volume of the shell, tube, cap and total
%clean up all previous data
clear
clc 

% Global Variables 
num_tubes = 6;
g = 9.81; 
% Case 1 Parameters
L1 = 3.000;   % Length of the shell in m
a1 = 0.400;   % Semi-major axis of the elliptical cross-section in m
b1 = 0.200;   % Semi-minor axis of the elliptical cross-section in m
w1 = 0.100;   % Width of the tube cross-section in m
d1 = 0.050;   % Depth of the tube cross-section in m
% Case 1 Parameters
L2 = 2.400;   % Length of the shell in m
a2 = 0.600;   % Semi-major axis of the elliptical cross-section in m
b2 = 0.400;   % Semi-minor axis of the elliptical cross-section in m
w2 = 0.150;   % Width of the tube cross-section in m
d2 = 0.075;   % Depth of the tube cross-section in m

% General notes: 
% All values in this code is in m unless stated. 
% This is because we have to everything in m

%-------------------------------------------------------------------------%
% Volume of shell 
% Performing the integration using Polar Coordinates
% Define sysmbolic variable
syms r theta y;  
% Rewriting z and x in the form of r and theta
% For case 1
% Define the function to integrate
f1_1 = r;
% Define the limits of integration
r_min1 = 0 ;
r_max1 = sqrt(1 ./ ((((cos(theta)).^2) ./ (a1.^2)) + (((sin(theta)).^2) ./ (b1.^2)))); 
y_min1 = 0;
y_max1 = L1;
theta_min1 = 0;
theta_max1 = 2*pi;
% Perform the numerical integration
% Integrating the first part with respect to z 
int_z_1 = int(f1_1, r, r_min1, r_max1);
int_x_1 = int(int_z_1, theta, theta_min1, theta_max1);  
F_part1_1 = int( int_x_1, y_min1, y_max1);
% For case 2
% Define the function to integrate
f1_2 = r;
% Define the limits of integration, converting x and y into pheta and r
r_min2 = 0 ;
r_max2 = 2*pi;
y_min2 = 0;
y_max2 = L2;
theta_min2 = 0;
theta_max2 = sqrt(1 ./ (((cos(theta).^2) ./ (a2.^2)) + ((sin(theta).^2) ./ (b2.^2))));
% Perform the numerical integration
% Integrating the first part with respect to z 
int_z_2 = int(f1_2, r , theta_min2, theta_max2);
int_x_2 = int(int_z_2, theta, r_min2, r_max2); 
F_part1_2 = int( int_x_2, y, y_min2, y_max2);
% Saving the values of the shells
% Convert symbolic to double
v_shell_1 = double(F_part1_1);
v_shell_2 = double(F_part1_2);
% Display the volume of the shell
disp(['The volume of the shell for case 1 is: ', num2str(v_shell_1), ' m^3']);
disp(['The volume of the shell for case 2 is: ', num2str(v_shell_2), ' m^3']);


%-------------------------------------------------------------------------%
% Volume of one shell a tube 

% Using triple integral to calculte the volume of tube 
% Define sysmbolic variable
syms x y z;

% For case 1
% Define the function to integrate
f3_1 = 1;
% We will calculate the integral based on previous question
init_dis1 = 0.150;  
spacing_bet_tubes_1 = 0.400; % Distance between the left most node and the next one (0.1 + 0.3)
volume_tubes_1 = 0; % Initializing the volume of the tubes as accumulation is supported in matlab
for i = 1:num_tubes 
    % Define the limits of integration
    x_min3_1 = - a1 ;
    x_max3_1 = a1;
    y_min3_1 = init_dis1 + spacing_bet_tubes_1.*(i-1);
    y_max3_1 = init_dis1 + spacing_bet_tubes_1.*(i-1) + w1;
    z_min3_1 = -d1./2;
    z_max3_1 = d1./2;
    % Perform the numerical integration
    % Integrating the first part with respect to z 
    int_z_3_1 = int(f3_1, z, z_min3_1, z_max3_1);
    int_x_3_1 = int(int_z_3_1, x_min3_1, x_max3_1); 
    F_part3_1 = int(int_x_3_1, y_min3_1, y_max3_1);
    volume_tubes_1 =  volume_tubes_1 + F_part3_1 ;
end

% For case 2
% Define the function to integrate
f3_2 = 1;
% We will calculate the integral based on previous question
init_dis2 = 0;  
spacing_bet_tubes_2 = 0.400; % Distance between the left most node and the next one
volume_tubes_2 = 0; 
% Initializing the volume of the tubes as accumulation is supported in matlab
for i = 1:num_tubes    
    % Define the limits of integration
    x_min3_2 = - a2 ;
    x_max3_2 = a2;
    y_min3_2 = init_dis2 + spacing_bet_tubes_2.*(i-1);
    y_max3_2 = init_dis2 + spacing_bet_tubes_2.*(i-1) + w2;
    z_min3_2 = -d2./2;
    z_max3_2 = d2./2;
    % Perform the numerical integration
    % Integrating the first part with respect to z 
    int_z_3_2 = int(f3_2, z, z_min3_2, z_max3_2);
    int_x_3_2 = int(int_z_3_2, x_min3_2, x_max3_2); 
    F_part3_2 = int(int_x_3_2, y_min3_2, y_max3_2);
    volume_tubes_2 =  volume_tubes_2 + F_part3_2 ;
end
% Saving the values of the tubes
% Convert symbolic to double
volume_tubes_1 = double(volume_tubes_1);
volume_tubes_2 = double(volume_tubes_2);
vol_tube_1 = volume_tubes_1 ./ 6; 
vol_tube_2 = volume_tubes_2 ./ 6; 
% Display the result for tubes 
disp(['The volume of one tube for case 1 is: ', num2str(vol_tube_1), ' m^3']);
disp(['The volume of one tube for case 2 is: ', num2str(vol_tube_2), ' m^3']);


%-------------------------------------------------------------------------%
% Volume of caps 
% Using Cartesian Coordinates 
num_tubes = 6 ; 
spacing_bet_tubes_1 = 0.4;
spacing_bet_tubes_2 = 0.35;
% For case 1
% Define the function to integrate
f4_1 = 1;
% We will calculate the integral based on previous question 1a
init_dis1 = 0.150; 
v_caps_1 = 0; % Initializing the volume of the tubes as accumulation is supported in matlab
for i = 1:num_tubes 
    % Define the limits of integration 
    x_min4_1 = sqrt((a1^2)*(1 - (z^2) / (b1^2))); 
    x_max4_1 = a1; 
    y_min4_1 = init_dis1 + spacing_bet_tubes_1.*(i-1);
    y_max4_1 = init_dis1 + spacing_bet_tubes_1.*(i-1) + w1;
    z_min4_1 = -d1./2;
    z_max4_1 = d1./2;
    % Perform the numerical integration for cap on the +ve x-axis
    % Integrating the first part with respect to z 
    int_x_4_1 = int(f4_1, x, x_min4_1, x_max4_1);
    int_z_4_1 = int(int_x_4_1, z, z_min4_1, z_max4_1); 
    F_part4_1 = int(int_z_4_1, y, y_min4_1, y_max4_1);
    v_caps_1 =  v_caps_1 + F_part4_1 ;
end
% Since there is caps on both side (1 on both the +ve and -ve x-axis), 
% Refer to diagrams in report (12 caps) 
v_caps_1 = v_caps_1 * 2; 

% For case 2
% Define the function to integrate % dx dz dy 
f4_2 = 1;
% We will calculate the integral based on previous question 1a
init_dis2 = 0.50;  
v_caps_2 = 0;
for i = 1:num_tubes 
    % Define the limits of integration 
    x_min4_2 = sqrt((a2^2)*(1 - (z^2) / (b2^2)));  
    x_max4_2 = a2;
    y_min4_2 = init_dis2 + spacing_bet_tubes_2.*(i-1);
    y_max4_2 = init_dis2 + spacing_bet_tubes_2.*(i-1) + w2;
    z_min4_2 = -d2./2;
    z_max4_2 = d2./2;
    % Perform the numerical integration for cap on the +ve x-axis
    % Integrating the first part with respect to z 
    int_x_4_2 = int(f4_2, x, x_min4_2, x_max4_2);
    int_z_4_2 = int(int_x_4_2, z, z_min4_2, z_max4_2); 
    F_part4_2 = int(int_z_4_2, y, y_min4_2, y_max4_2);
    v_caps_2 =  v_caps_2 + F_part4_2 ;
end 
% Since there is caps on both side (1 on both the +ve and -ve x-axis), 
% Refer to diagrams in report (12 caps) 
v_caps_2 = v_caps_2 * 2; 

% Saving the values of the tubes
% Convert symbolic to double
v_caps_1 = double(v_caps_1);
v_caps_2 = double(v_caps_2);
v_1caps_1 = v_caps_1 ./ 12 ; 
v_1caps_2 = v_caps_2 ./ 12 ; 

% Display the volume of the shell
disp(['The volume of one cap for case 1 is: ', num2str(v_1caps_1), ' m^3']);
disp(['The volume of one cap for case 2 is: ', num2str(v_1caps_2), ' m^3']);


%-------------------------------------------------------------------------%
% Total volume of oil 
% Based on the formula given in the project description
v_total_1 = v_shell_1 + v_caps_1 - 6 * vol_tube_1;
v_total_2 = v_shell_2 + v_caps_2 - 6 * vol_tube_2;
disp(['The total volume of oil in heat exchanger for case 1 is: ', num2str(v_total_1), ' m^3']);
disp(['The total volume of oil in heat exchanger for case 2 is: ', num2str(v_total_2), ' m^3']);

%-------------------------------------------------------------------------%
% Heat Exchanger Center of Mass

% General equation of the Density of Oil in the Shell
% Defining this as a general equation as it is used several time
% The '(@theta)' is used to ensure that correct sysmbolic variable is used
% and that the theta is treated as a variable. 
% This is important as I was having weird outputs when the (@theta) was not
% included. This is mainly because the formula is used multiple times. 
formula_density_oil = @(theta) 420 + 80 * (cos(theta - pi/2) + 5 * sin(theta + 5*pi/2));

% Mass of the Shell in kg 
% Performing the integration using Polar Coordinates
% Define sysmbolic variable
syms r theta y;  
% Rewriting z and x in the form of r and theta
% For case 1
% Define the function to integrate
f5_1 = formula_density_oil(theta);
% Define the limits of integration
r_min5_1 = 0 ;
r_max5_1 = sqrt(1 ./ ((((cos(theta)).^2) ./ (a1.^2)) + (((sin(theta)).^2) ./ (b1.^2)))); 
y_min5_1 = 0;
y_max5_1 = L1;
theta_min5_1 = 0;
theta_max5_1 = 2*pi;
% Perform the numerical integration
% Integrating the first part with respect to z 
int_z_5_1 = int(f5_1*r, r, r_min5_1, r_max5_1);
int_x_5_1 = int(int_z_5_1, theta, theta_min5_1, theta_max5_1);  
F_part5_1 = int(int_x_5_1, y, y_min5_1, y_max5_1);
% Mass of Oil in Shell for case 1
mass_oil_1 = double(F_part5_1);
% For case 2
% Define the function to integrate
f5_2 = formula_density_oil(theta);
% Define the limits of integration
r_min5_2 = 0 ;
r_max5_2 = sqrt(1 ./ ((((cos(theta)).^2) ./ (a2.^2)) + (((sin(theta)).^2) ./ (b2.^2)))); 
y_min5_2 = 0;
y_max5_2 = L2;
theta_min5_2 = 0;
theta_max5_2 = 2*pi;
% Perform the numerical integration
% Integrating the first part with respect to z 
int_z_5_2 = int(f5_2*r, r, r_min5_2, r_max5_2);
int_x_5_2 = int(int_z_5_2, theta, theta_min5_2, theta_max5_2);  
F_part5_2 = int(int_x_5_2, y, y_min5_2, y_max5_2);
% Mass of Oil in Shell for case 2
mass_oil_2 = double(F_part5_2);

% Display the mass of shell in kg
disp(['The mass of shell in kg for case 1 is: ', num2str(mass_oil_1), ' kg']);
disp(['The mass of shell in kg for case 2 is: ', num2str(mass_oil_2), ' kg']);


%-------------------------------------------------------------------------%
% Center of mass in x and z
% Mass of the Shell in kg 
% Performing the integration using Polar Coordinates
% Define sysmbolic variable
syms r theta y; 
% General equations of the center of mass in the shell
x_norm = r * cos(theta);
z_norm = r * sin(theta);
m1 = mass_oil_1; 
m2 = mass_oil_2; 

% For case 1
% Define the function to integrate
f6_1 = formula_density_oil(theta);
% Define the limits of integration
r_min6_1 = 0 ;
r_max6_1 = sqrt(1 ./ ((((cos(theta)).^2) ./ (a1.^2)) + (((sin(theta)).^2) ./ (b1.^2)))); 
y_min6_1 = 0;
y_max6_1 = L1;
theta_min6_1 = 0;
theta_max6_1 = 2*pi;
% Perform the numerical integration
% Integrating to find the center of x; moment of inertia 
int_z_6_1_1 = int(f6_1*x_norm*r, r, r_min6_1, r_max6_1);
int_x_6_1_1 = int(int_z_6_1_1, theta, theta_min6_1, theta_max6_1);  
F_part6_1_x_1 = int(int_x_6_1_1, y, y_min6_1, y_max6_1);
% Integrating to find the center of z; moment of inertia 
int_z_6_1_2 = int(f6_1*z_norm*r, r, r_min6_1, r_max6_1);
int_x_6_1_2 = int(int_z_6_1_2, theta, theta_min6_1, theta_max6_1);  
F_part6_1_z_1 = int(int_x_6_1_2, y, y_min6_1, y_max6_1);
% Center of mass in the x and z axis for case 1
center_mass_x_1 = double(F_part6_1_x_1./m1);
center_mass_z_1 = double(F_part6_1_z_1./m1);

% For case 2
% Define the function to integrate
f6_2 = formula_density_oil;
% Define the limits of integration
r_min6_2 = 0 ;
r_max6_2 = sqrt(1 ./ ((((cos(theta)).^2) ./ (a2.^2)) + (((sin(theta)).^2) ./ (b2.^2)))); 
y_min6_2 = 0;
y_max6_2 = L2;
theta_min6_2 = 0;
theta_max6_2 = 2*pi;
% Perform the numerical integration
% Integrating to find the center of x 
int_z_6_2_1 = int(f6_2*x_norm*r, r, r_min6_2, r_max6_2);
int_x_6_2_1 = int(int_z_6_2_1, theta, theta_min6_2, theta_max6_2);  
F_part6_2_x_1 = int(int_x_6_2_1, y, y_min6_2, y_max6_2);
% Integrating to find the center of z 
int_z_6_2_2 = int(f6_2*z_norm*r, r, r_min6_2, r_max6_2);
int_x_6_2_2 = int(int_z_6_2_2, theta, theta_min6_2, theta_max6_2);  
F_part6_2_z_1 = int(int_x_6_2_2, y, y_min6_2, y_max6_2);
% Center of mass in the x and z axis for case 2
center_mass_x_2 = double(F_part6_2_x_1./m2);
center_mass_z_2 = double(F_part6_2_z_1./m2);


% Display the center of mass in the x and z axis
disp(['The center of mass in the x-axis for case 1 is: ', num2str(center_mass_x_1), ' m']);
disp(['The center of mass in the z-axis for case 1 is: ', num2str(center_mass_z_1), ' m']);
disp(['The center of mass in the x-axis for case 2 is: ', num2str(center_mass_x_2), ' m']);
disp(['The center of mass in the z-axis for case 2 is: ', num2str(center_mass_z_2), ' m']);


%-------------------------------------------------------------------------%
% Torque along the x-axis from point A

% Case 1 
Torque_1 = m1*g*(center_mass_x_1+a1); 
% Case 2 
Torque_2 = m2*g*(center_mass_x_2+a2); 

% Display the torque of the two cases
disp(['The torque along the x-axis for case 1 is: ', num2str(Torque_1), ' Nm']);
disp(['The torque along the x-axis for case 2 is: ', num2str(Torque_2), ' Nm']);