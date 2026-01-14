%% trajectory_spiral.m
% Definizione della traiettoria desiderata per il quadcopter

% Parametri spirale
R = 1.0;            % raggio [m]
omega_traj = 0.5;   % velocità angolare [rad/s]
z0 = 0.0;           % quota iniziale
zf = 1.5;           % quota finale

% Funzione di velocità verticale
Tend = 10;          % durata simulazione (da sincronizzare)
vz = (zf - z0)/Tend;

% Funzione posizione desiderata
xi_d = @(tt) [ R*cos(omega_traj*tt);
               R*sin(omega_traj*tt);
               z0 + vz*tt ];

% Funzione velocità
xi_d_dot = @(tt) [ -R*omega_traj*sin(omega_traj*tt);
                    R*omega_traj*cos(omega_traj*tt);
                    vz ];

% Funzione accelerazione
xi_d_ddot = @(tt) [ -R*omega_traj^2*cos(omega_traj*tt);
                    -R*omega_traj^2*sin(omega_traj*tt);
                     0 ];