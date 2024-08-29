% Part 1 a)
% Using MATLAB to plot the tank and the 6 tubes for the heat exchanger
% Plotted the case 1 using the Table 1

% Clearing previous data for safety 
clear; 
clc; 

% Define parameters for the main cylinder
% Case 1 Parameters
L1 = 3.000;   % Length of the shell in m
a1 = 0.400;   % Semi-major axis of the elliptical cross-section in m
b1 = 0.200;   % Semi-minor axis of the elliptical cross-section in m
w1 = 0.100;   % Width of the tube cross-section in m
d1 = 0.050;   % Depth of the tube cross-section in m

% A bit redundant but wanted to be able to use the code faster for
% different values of L, a, b, w and d 
L = L1;  
a = a1;   
b = b1;   
w = w1;   
d = d1;    

% Creating the shape of the shell; defining the bounds of theta and y 
theta = linspace(0, 2*pi, 100); 
y = linspace(0, L, 100);
[Theta, Y] = meshgrid(theta, y);

% Parametrization
X_shell = a * cos(Theta);
Z_shell = b * sin(Theta);
Y_shell = Y;

% Plot the main cylinder
figure;
surf(X_shell, Y_shell, Z_shell, 'FaceAlpha', 0.5, 'EdgeColor', 'none');
hold on; % Hold on to be able to do the tubes later

% Create the Rectangular tubes, got help from ChatGPT and the TA
num_tubes = 6; 
tube_offset_x = w/2; 
tube_offset_y = w ; 
tube_positions = linspace(-L/2, L/2, num_tubes); % Adjusted for X-axis
spacing_tubes = 0.300 + w;  
% 0.300 is calculated by looking for the nearest factor of 7 (num of tubes + 1)
% which is 2.100 here and 2.100/7 = 0.300

% Define positions for the small tubes in the y-z plane
% The help of a TA and ChatGPT was used for this section
% A for loop is used as there is no change in my my x and z values but 
for i = 1:num_tubes
    % Here we used the idea of plotting the different points on the
    % 3 different axis
    x_tube = [-((2*a)/2), ((2*a)/2), ((2*a)/2), -((2*a)/2), -((2*a)/2), ((2*a)/2), ((2*a)/2), -((2*a)/2)] - tube_offset_x;
    y_tube = [-w/2, -w/2, -w/2, -w/2, w/2, w/2, w/2, w/2] + tube_offset_y + (i)*spacing_tubes;
    z_tube = [-d/2, -d/2, d/2, d/2, -d/2, -d/2, d/2, d/2] + d/2 - 0.275;
    
    % Then we put the tube to the correct position
    y_tube_center = tube_offset_x;
    z_tube_center = a/2 + w/2;
    tube = [x_tube' + y_tube_center, y_tube', z_tube'+ z_tube_center];
    
    % Define the faces
    faces = [1 2 3 4; 5 6 7 8; 1 2 6 5; 2 3 7 6; 3 4 8 7; 4 1 5 8]; %Front; back faces

 % Plot the tube
    patch('Vertices', tube, 'Faces', faces, 'FaceColor', 'red', 'FaceAlpha', 0.5);
    
end

% Adjust plot settings
xlabel('X-axis');
ylabel('Y-axis');
zlabel('Z-axis');
title('Case 1: Tank with 6 Equidistant Tubes')
% title('Case 2: Tank with 6 Equidistant Tubes')

axis equal;
% view([1 0 0]); % This sets the view to the z-y plane
% Was used to explain how the calibrating was done




