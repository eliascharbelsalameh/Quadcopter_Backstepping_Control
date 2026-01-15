%% trajectory_lemniscate.m

% Parametri lemniscata
A = 1.0;          % Amplitude 
omega_traj = 0.5; % angular speed [rad/s]
z0 = 0.0;         % initial altitude
zf = 1.5;         % final altitude
Tend = 10;        % duration
vz = (zf - z0)/Tend;  % vertical velocity

% Lemnicate function
xi_d = @(tt) [ A * cos(omega_traj*tt) ./ (1 + sin(omega_traj*tt).^2);  % x
               A * sin(omega_traj*tt).*cos(omega_traj*tt) ./ (1 + sin(omega_traj*tt).^2);  % y
               z0 + vz*tt ];  % z

% I derivate
xi_d_dot = @(tt) [ -A*omega_traj*sin(omega_traj*tt)./(1 + sin(omega_traj*tt).^2) ...
                    - A*cos(omega_traj*tt).* (2*sin(omega_traj*tt)*omega_traj.*cos(omega_traj*tt)) ./ (1 + sin(omega_traj*tt).^2).^2;
                    
                    A*omega_traj*(cos(omega_traj*tt).^2 - sin(omega_traj*tt).^2) ./ (1 + sin(omega_traj*tt).^2) ...
                    - A*sin(omega_traj*tt).*cos(omega_traj*tt) .* (2*sin(omega_traj*tt)*omega_traj.*cos(omega_traj*tt)) ./ (1 + sin(omega_traj*tt).^2).^2;
                    
                    vz ];

% II derivate
xi_d_ddot = @(tt) [ ...
    -A*omega_traj^2 * cos(omega_traj*tt)./(1 + sin(omega_traj*tt).^2) ...
    + 4*A*omega_traj^2*sin(omega_traj*tt).^2.*cos(omega_traj*tt).^2./(1 + sin(omega_traj*tt).^2).^3 ...
    - 2*A*omega_traj^2*sin(omega_traj*tt).*cos(omega_traj*tt).^3./(1 + sin(omega_traj*tt).^2).^3;
    
    -2*A*omega_traj^2*sin(omega_traj*tt).*cos(omega_traj*tt)./(1 + sin(omega_traj*tt).^2) ...
    + 2*A*omega_traj^2*(cos(omega_traj*tt).^2 - sin(omega_traj*tt).^2).*sin(omega_traj*tt).*cos(omega_traj*tt)./(1 + sin(omega_traj*tt).^2).^2 ...
    - 4*A*omega_traj^2*sin(omega_traj*tt).^2.*cos(omega_traj*tt).^2./(1 + sin(omega_traj*tt).^2).^3;
    
    0 ];
