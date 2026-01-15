%% trajectory_spiral.m
% Definition of the desired trajectory for the quadcopter (Conical Spiral)

% Spiral Parameters
R_start    = 0.5;          % Initial radius [m]
R_end      = 2.5;          % Final radius [m]
omega_traj = 0.8;          % Angular velocity [rad/s]
z0         = 0.0;          % Initial altitude
zf         = 2.0;          % Final altitude

% Simulation duration 
Tend       = 40;           

% Growth rates
vz = (zf - z0) / Tend;          % Vertical velocity (constant)
vr = (R_end - R_start) / Tend;  % Radial expansion velocity (constant)

% 1. Desired Position Function
% r(t) = R_start + vr * t
xi_d = @(tt) [ (R_start + vr*tt) * cos(omega_traj*tt);
               (R_start + vr*tt) * sin(omega_traj*tt);
               z0 + vz*tt ];

% 2. Desired Velocity Function
% Uses Product Rule: d/dt(r*cos) = r'*cos + r*(cos)'
xi_d_dot = @(tt) [ vr*cos(omega_traj*tt) - (R_start + vr*tt)*omega_traj*sin(omega_traj*tt);
                   vr*sin(omega_traj*tt) + (R_start + vr*tt)*omega_traj*cos(omega_traj*tt);
                   vz ];

% 3. Desired Acceleration Function
% Uses Product Rule again on velocity terms
% Contains Centripetal term (R*w^2) and Coriolis terms (2*vr*w)
xi_d_ddot = @(tt) [ -2*vr*omega_traj*sin(omega_traj*tt) - (R_start + vr*tt)*omega_traj^2*cos(omega_traj*tt);
                     2*vr*omega_traj*cos(omega_traj*tt) - (R_start + vr*tt)*omega_traj^2*sin(omega_traj*tt);
                     0 ];