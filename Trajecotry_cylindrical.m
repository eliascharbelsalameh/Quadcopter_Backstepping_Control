%% trajectory_cylindrical.m
% Traiettoria desiderata: cilindrica
% Definizione posizione, velocità e accelerazione

% Parametri cilindro
R = 1.0;          % raggio del cilindro [m]
omega_traj = 0.5; % velocità angolare attorno al cilindro [rad/s]
z0 = 0.0;         % quota iniziale
zf = 1.5;         % quota finale
Tend = 10;        % durata simulazione (da sincronizzare con lo script principale)
vz = (zf - z0)/Tend;  % velocità verticale costante

% Funzione posizione desiderata (cilindro verticale)
xi_d = @(tt) [ R * cos(omega_traj*tt);   % x
               R * sin(omega_traj*tt);   % y
               z0 + vz*tt ];             % z

% Velocità desiderata (derivata prima)
xi_d_dot = @(tt) [ -R * omega_traj * sin(omega_traj*tt);  % ẋ
                    R * omega_traj * cos(omega_traj*tt);  % ẏ
                    vz ];                                  % ż

% Accelerazione desiderata (derivata seconda)
xi_d_ddot = @(tt) [ -R * omega_traj^2 * cos(omega_traj*tt);  % ẍ
                    -R * omega_traj^2 * sin(omega_traj*tt);  % ÿ
                    0 ];                                      % z̈
