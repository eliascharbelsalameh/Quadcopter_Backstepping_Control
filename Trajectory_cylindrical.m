%% trajectory_cylindrical.m


% Parametri cilindro
R = 1.0;          % radius [m]
omega_traj = 0.5; % angular velocity [rad/s]
z0 = 0.0;         % initial altitude
zf = 1.5;         % final altitude
Tend = 10;        % duration
vz = (zf - z0)/Tend;  % vertical velocity

% Desiderd function
xi_d = @(tt) [ R * cos(omega_traj*tt);   % x
               R * sin(omega_traj*tt);   % y
               z0 + vz*tt ];             % z

% I derivate
xi_d_dot = @(tt) [ -R * omega_traj * sin(omega_traj*tt);  
                    R * omega_traj * cos(omega_traj*tt);  
                    vz ];                               ̇

% II derivate
xi_d_ddot = @(tt) [ -R * omega_traj^2 * cos(omega_traj*tt);  
                    -R * omega_traj^2 * sin(omega_traj*tt); 
                    0 ];                                      
