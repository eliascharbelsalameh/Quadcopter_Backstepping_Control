%% point_target.m

% Target point
target_point = [10.0; 10.0; 10];  % x, y, z

% desired function and its derivatives
xi_d = @(tt) target_point;       
xi_d_dot = @(tt) [0;0;0];        
xi_d_ddot = @(tt) [0;0;0];       
