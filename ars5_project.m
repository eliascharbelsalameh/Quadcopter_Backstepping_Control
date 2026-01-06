%% ARS5 Project
% authors: 
%           Stefano MUSSO PIANTELLI
%           Elias Charbel SALAMEH
% issued on: 06/01/2026

clear; 
clc; 
close all;

%% Setup
g = 9.81;
Tend = 40;          
dt   = 5e-3;        % Slightly smaller dt for filter accuracy
t    = 0:dt:Tend;   
N    = numel(t);    

%% Controller Tuning
k1 = 2.0; % Controller Gains
k2 = 5.0; 
k3 = 5.0; 
k4 = 5.0; 

% DSC Filter Time Constant (tau)
% Determines how fast the filter tracks the raw signal.
% Too small = noisy derivative. Too large = delay/instability.
tau = 0.02; 

%% Reference Trajectory Definition 
x1_star = @(tt) exp(-0.2*tt).*sin(3*tt) + 0.5*sin(0.5*tt);

x1_star_dot = @(tt) exp(-0.2*tt).*(3*cos(3*tt) - 0.2*sin(3*tt)) ...
    + 0.25*cos(0.5*tt);

x1_star_ddot = @(tt) exp(-0.2*tt).*(-8.96*sin(3*tt) - 1.2*cos(3*tt)) ...
    - 0.125*sin(0.5*tt);

%% Memory Allocation & Initialization
x = zeros(4, N);
u = zeros(1, N);
e = zeros(4, N);

% phi_f  : Filtered Virtual Controls (Actual values used in error calc)
% phi_dot: Smooth derivatives from the filter
phi_f   = zeros(4, N);
phi_dot = zeros(4, N);

% Initial Conditions
x(:,1) = [0; 0; 0; 0];

% Initialize filters with steady state assumptions to avoid startup "jumps"
% At t=0, x1_star_dot is approx 3.25. 
% We initialize phi_f(2) to this so it doesn't start at 0 and jump.
phi_f(2,1) = x1_star_dot(0); 

%% Main Simulation Loop
fprintf('Starting DSC simulation...\n');

for k = 1:N-1
    tk = t(k);
    
    % --- Step 1: Position Error ---
    x1d     = x1_star(tk);
    x1d_dot = x1_star_dot(tk);
    
    e(1,k) = x(1,k) - x1d;
    
    % Raw Virtual Control 2 (Target)
    phi2_cmd = -k1 * e(1,k) + x1d_dot;
    
    % DSC Filter for phi2
    % Dynamics: tau * phi_dot + phi = phi_cmd
    % Therefore: phi_dot = (phi_cmd - phi_current) / tau
    phi_dot(2,k) = (phi2_cmd - phi_f(2,k)) / tau;
    
    % Update Filter State (Euler Integration)
    phi_f(2,k+1) = phi_f(2,k) + dt * phi_dot(2,k);
    
    
    % --- Step 2: Velocity Error ---
    % Note: We use the FILTERED value phi_f(2) for the error, not the raw command
    e(2,k) = x(2,k) - phi_f(2,k);
    
    % Raw Virtual Control 3 (Target)
    % Inverting x2_dot = -g * tan(x3)
    term_inside = (k2 * e(2,k) - phi_dot(2,k)) / g;
    phi3_cmd    = atan(term_inside);
    
    % DSC Filter for phi3
    phi_dot(3,k) = (phi3_cmd - phi_f(3,k)) / tau;
    
    % Update Filter State
    phi_f(3,k+1) = phi_f(3,k) + dt * phi_dot(3,k);
    
    
    % --- Step 3: Angle Error ---
    e(3,k) = x(3,k) - phi_f(3,k);
    
    % Raw Virtual Control 4 (Target)
    phi4_cmd = -k3 * e(3,k) + phi_dot(3,k);
    
    % DSC Filter for phi4
    phi_dot(4,k) = (phi4_cmd - phi_f(4,k)) / tau;
    
    % Update Filter State
    phi_f(4,k+1) = phi_f(4,k) + dt * phi_dot(4,k);
    
    
    % --- Step 4: Angle Rate Error ---
    e(4,k) = x(4,k) - phi_f(4,k);
    
    % Final Control Law (No filter needed here, this is the actual input)
    u(k) = -k4 * e(4,k) + phi_dot(4,k);
    
    
    % --- Step 5: System Dynamics ---
    % Safety Saturation
    if abs(x(3,k)) > 1.56 
       x(3,k) = sign(x(3,k)) * 1.56; 
    end

    x1_dot = x(2,k);
    x2_dot = -g * tan(x(3,k));
    x3_dot = x(4,k);
    x4_dot = u(k);
    
    % Update Plant State
    x(1,k+1) = x(1,k) + dt * x1_dot;
    x(2,k+1) = x(2,k) + dt * x2_dot;
    x(3,k+1) = x(3,k) + dt * x3_dot;
    x(4,k+1) = x(4,k) + dt * x4_dot;
end

% Final values for plotting
u(end) = u(end-1);
e(:,end) = e(:,end-1);
fprintf('Simulation complete.\n');

%% Visualization
figure('Color','w','Name','DSC Backstepping Results', 'Position', [100 100 1000 700]);

% Trajectory
subplot(3,2,1:2);
plot(t, x(1,:), 'b', 'LineWidth', 1.5); hold on;
plot(t, arrayfun(x1_star, t), 'r--', 'LineWidth', 1.5);
ylabel('Position x_1'); title('Trajectory Tracking (DSC)');
legend('Actual','Reference'); grid on;

% Angle
subplot(3,2,3);
plot(t, x(3,:), 'Color', [0.8 0.4 0], 'LineWidth', 1.2);
ylabel('Angle x_3 [rad]'); title('System Angle');
yline(pi/2,'y:'); yline(-pi/2,'y:'); grid on;

% Control Input (Notice how much smoother this is compared to finite diff!)
subplot(3,2,4);
plot(t, u, 'c', 'LineWidth', 1.0);
ylabel('u'); title('Control Input (Smoothed)');
grid on;

% Errors
subplot(3,2,5:6);
plot(t, e(1,:), 'LineWidth', 1.2); hold on;
plot(t, e(2,:), 'LineWidth', 1.2);
plot(t, e(3,:), 'LineWidth', 1.2);
plot(t, e(4,:), 'LineWidth', 1.2);
ylabel('Error'); xlabel('Time [s]');
legend('e1','e2','e3','e4'); title('Tracking Errors');
grid on;