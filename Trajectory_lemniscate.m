%% trajectory_lemniscate.m
% Traiettoria desiderata: lemniscata (∞)
% Definizione posizione, velocità e accelerazione

% Parametri lemniscata
A = 1.0;         % ampiezza (semi-asse principale)
omega_traj = 0.5; % velocità angolare [rad/s]
z0 = 0.0;         % quota iniziale
zf = 1.5;         % quota finale
Tend = 10;        % durata simulazione (da sincronizzare con lo script principale)
vz = (zf - z0)/Tend;  % velocità verticale costante

% Funzione posizione desiderata (lemniscata di Bernoulli)
xi_d = @(tt) [ A * cos(omega_traj*tt) ./ (1 + sin(omega_traj*tt).^2);  % x
               A * sin(omega_traj*tt).*cos(omega_traj*tt) ./ (1 + sin(omega_traj*tt).^2);  % y
               z0 + vz*tt ];  % z

% Velocità desiderata (derivata simbolica)
xi_d_dot = @(tt) [ -A*omega_traj*sin(omega_traj*tt)./(1 + sin(omega_traj*tt).^2) ...
                    - A*cos(omega_traj*tt).* (2*sin(omega_traj*tt)*omega_traj.*cos(omega_traj*tt)) ./ (1 + sin(omega_traj*tt).^2).^2;
                    
                    A*omega_traj*(cos(omega_traj*tt).^2 - sin(omega_traj*tt).^2) ./ (1 + sin(omega_traj*tt).^2) ...
                    - A*sin(omega_traj*tt).*cos(omega_traj*tt) .* (2*sin(omega_traj*tt)*omega_traj.*cos(omega_traj*tt)) ./ (1 + sin(omega_traj*tt).^2).^2;
                    
                    vz ];

% Accelerazione desiderata (derivata seconda)
xi_d_ddot = @(tt) [ ...
    -A*omega_traj^2 * cos(omega_traj*tt)./(1 + sin(omega_traj*tt).^2) ...
    + 4*A*omega_traj^2*sin(omega_traj*tt).^2.*cos(omega_traj*tt).^2./(1 + sin(omega_traj*tt).^2).^3 ...
    - 2*A*omega_traj^2*sin(omega_traj*tt).*cos(omega_traj*tt).^3./(1 + sin(omega_traj*tt).^2).^3;
    
    -2*A*omega_traj^2*sin(omega_traj*tt).*cos(omega_traj*tt)./(1 + sin(omega_traj*tt).^2) ...
    + 2*A*omega_traj^2*(cos(omega_traj*tt).^2 - sin(omega_traj*tt).^2).*sin(omega_traj*tt).*cos(omega_traj*tt)./(1 + sin(omega_traj*tt).^2).^2 ...
    - 4*A*omega_traj^2*sin(omega_traj*tt).^2.*cos(omega_traj*tt).^2./(1 + sin(omega_traj*tt).^2).^3;
    
    0 ];
