%% ARS5 Project
% Title:    Quadcopter Backstepping Controller
% 
% Authors:  Stefano MUSSO PIANTELLI, 
%           Elias Charbel SALAMEH
% 
% Date:     12/01/2026

clear; 
clc; 
close all;

%% Parameters
% Physical parameters
g       = 9.81;
m       = 0.6;
Ixx     = 0.0014; Iyy = 0.0165; Izz = 0.0152;
I       = diag([Ixx Iyy Izz]);
F_ext   = [-0.1*randn(1,1);-0.1*randn(1,1);-g];

% Simulation parameters
Tend    = 20;
dt      = 1e-3;        % smaller dt for better attitude dynamics resolution
t       = 0:dt:Tend;
N       = numel(t);

%% Desired trajectory – Ascending spiral
% run('Trajectory_spiral.m');
% run('Trajectory_lemniscate.m');
run('Trajectory_cylindrical.m');
% run('Target_point.m');

%% State vector: X = [xi(3); v(3); q(4); Omega(3)]
% quaternion q = [qw; qx; qy; qz], unit norm
X       = zeros(13, N);

% Initial conditions:
xi0     = [10; -5; 7];                % initial position
%xi0     = xi_d(0);                  % instead of [0;0;0]
v0      = [2; 0; 2];                % initial velocity
q0      = [1; 0; 0; 0];             % initial attitude (no rotation)
Omega0  = [0; 0; 0];                % initial angular velocity
X(:,1)  = [xi0; v0; q0; Omega0];

%% Gains
k1  = 1;
k2  = 1;
k3  = 1;
k4  = 1;
a   = 1;
d   = 0;

%% New system
e1 = zeros(3, N);
e1_dot = zeros(3, N);
e2 = zeros(3, N);
e2_dot = zeros(3, N);
q_e = zeros(4, N);
alpha_2 = zeros(3, N);

xi_v_dot = zeros(3,N);
xi_v_ddot = zeros(3,N);
omega_v = zeros(3,N);
omega_v_dot = zeros(3,N);
q_e_dot = zeros(4,N);
alpha_2_dot = zeros(3,N);
Fu = zeros(3,N);
psi_d = zeros(1,N);

%% sampled loop
for k = 1:N-1
    tk = t(k);

    xi_k = X(1:3,k);
    xi_dot_k = X(4:6,k);
    q_k = X(7:10,k);
    omega_k = X(11:13,k);

    % target
    xi_d_val      = xi_d(tk);
    xi_d_dot_val  = xi_d_dot(tk);
    xi_d_ddot_val = xi_d_ddot(tk);




    e1(:,k) = -(a*(xi_k - xi_d_val)-d); %e1
    xi_v_dot(:,k) = (k1*e1(:,k))/a + xi_d_dot_val; %equation 23
    e2(:,k) = xi_dot_k - xi_v_dot(:,k); %e2
    e1_dot(:,k) = -a*e2(:,k)-k1*e1(:,k); %e1dot
    % F_des: control input
    vec1 = m * ( (a * e1(:,k)) - F_ext + ...
                 ((k1 * e1_dot(:,k)) / a + xi_d_ddot_val) - ...
                 (k2 * e2(:,k)) ); %equazione 30, qFq*


    q_row = q_k.'; %q as a row 
    q_left = quatinv(q_row); 
    % q_conj = quatconj(q_row);
    % q_right = quatinv(q_conj);
    vec_quat = [0, vec1'];

    %comput F_u
    F_des_quat = quatmultiply(quatmultiply(q_left, vec_quat), q_row);
    F_des = F_des_quat(2:4)';

    Fu(:,k) = F_des;
    F_des_norm = norm(F_des);

    % e2_dot
    vec2 = Fu(:,k)/m;
    q_conj2 = quatconj(q_row);
    vec_quat2 = [0, vec2'];
    e2_dot_quat = quatmultiply(quatmultiply(q_row, vec_quat2), q_conj2) + F_ext - (k1*e1_dot(:,k)/a + xi_d_ddot_val);
    e2_dot(:,k) = e2_dot_quat(2:4);
    % attitude extraction
    if F_des_norm < 1e-6
        % Fallback: Zero thrust and identity orientation (or keep previous)
        Fu_des_norm = 0;
        q_d = [1, 0, 0, 0]; % Identity quaternion [w x y z]
        warning('Fu magnitude is near zero. Thrust set to 0, attitude set to identity.');
        % tau = 0;
    end
    %z_d = F_des / Fu(:,k); %%%%your're doing a divison between vector, no sense
    z_d = F_des / F_des_norm;
    l_vec = [-sin(psi_d(:,k)); cos(psi_d(:,k)); 0];

    l_cross_zd = cross(l_vec, z_d);
    norm_l_cross_zd = norm(l_cross_zd);
    if norm_l_cross_zd < 1e-6
        if abs(z_d(3)) < 0.999
            temp_v = [0; 0; 1];
            x_d_temp = cross(temp_v, z_d);
            x_d = x_d_temp / norm(x_d_temp);
        else
            x_d = [1; 0; 0];
        end
        warning('Singularity detected in attitude extraction (l parallel to z_d). Using fallback x_d.');
    else
        x_d = l_cross_zd / norm_l_cross_zd;
    end    
    y_d = cross(z_d, x_d);
    Rd = [x_d, y_d, z_d];
    q_d = rotm2quat(Rd);
    % quaternion error
    q_k   = q_k / norm(q_k);
    q_d = q_d / norm(q_d);
    q_d_conj = quatconj(q_d);
    q_bar_e  = quatmultiply(q_d_conj, q_k'); %aggiunto un trasposto
    % q_e from q_bar_e
    q_bar_e0 = q_bar_e(1); q_bar_ev = q_bar_e(2:4).';
    q_e(:,k) = [1 - abs(q_bar_e0); q_bar_ev];
    % M : quaternion error dynamics computation
    M = quatM(q_bar_e);
    omega_v(:,k) = -2*k3*M.'*q_e(:,k);
    alpha_2(:,k) = omega_k - omega_v(:,k);
    q_e_dot(:,k) = 0.5 * M * (alpha_2(:, k) + omega_v(:,k)); %dimension
    M_dot = quatM(q_e_dot(:,k));
    omega_v_dot = -2*k3*(M_dot.' * q_e(:,k) + M.' * q_e_dot(:,k));
    % tau
    tau = cross(omega_k,I*omega_k)+I*omega_v_dot-k4*alpha_2(:, k)-0.5*M.'*q_e;
    alpha_2_dot(:,k) = I\(-cross(omega_k,I*omega_k)+tau(:, k))-omega_v_dot;

    % initial system
    xi_dot = xi_dot_k;
    
    % q_k = quatmultiply(quatinv(q_conj),q_e(:,k));
    % q_k = q_k / norm(q_k);
    %xi_ddot = quatmultiply(quatmultiply(q_row, vec_quat2), q_conj2) + F_ext;
    vec_rotated = quatmultiply(quatmultiply(q_row, vec_quat2), q_conj2);  % 1x4
    xi_ddot = vec_rotated(2:4).' + F_ext;  % 3x1
    q_dot = 0.5*quatmultiply(q_k',[0 (omega_k.')]);

    omega_dot = I\(tau-cross(omega_k,I*omega_k));

    % update
    X(1:3,k+1) = xi_k       + dt * xi_dot;
    X(4:6,k+1) = xi_dot_k   + dt * xi_ddot;
    X(7:10,k+1) = q_k       + dt * q_dot.';
    X(11:13,k+1) = omega_k  + dt * omega_dot(:, k);
end


%% Calcolo traiettoria desiderata completa
xi_d_all = zeros(3, N);  % 3 righe: x, y, z
for k = 1:N
    xi_d_all(:,k) = xi_d(t(k));
end

%% Calcolo traiettoria desiderata completa
xi_d_all = zeros(3, N);  % 3 righe: x, y, z
for k = 1:N
    xi_d_all(:,k) = xi_d(t(k));  % funzione xi_d restituisce colonna 3x1
end

%% Plotting dei Risultati

% Estrazione posizioni simulate dallo stato X
xi_sim = X(1:3, :);

% 1. Visualizzazione 3D della Traiettoria
figure('Name', 'Traiettoria Quadricottero 3D', 'Color', 'w');
plot3(xi_d_all(1,:), xi_d_all(2,:), xi_d_all(3,:), 'g--', 'LineWidth', 1.5); hold on;
plot3(xi_sim(1,:), xi_sim(2,:), xi_sim(3,:), 'c', 'LineWidth', 1.5);
grid on; axis equal;
xlabel('X [m]'); ylabel('Y [m]'); zlabel('Z [m]');
legend('Desiderata (Spiral)', 'Simulata (Backstepping)');
title('Inseguimento Traiettoria 3D');
view(45, 30);

% 2. Analisi temporale per ogni asse (x, y, z)
figure('Name', 'Analisi Assi x, y, z', 'Color', 'w');

subplot(3,1,1);
plot(t, xi_d_all(1,:), 'g--', t, xi_sim(1,:), 'c');
ylabel('x [m]'); grid on;
legend('Desiderata', 'Simulata');
title('Inseguimento Posizione sui singoli assi');

subplot(3,1,2);
plot(t, xi_d_all(2,:), 'g--', t, xi_sim(2,:), 'c');
ylabel('y [m]'); grid on;

subplot(3,1,3);
plot(t, xi_d_all(3,:), 'g--', t, xi_sim(3,:), 'c');
ylabel('z [m]'); xlabel('Tempo [s]'); grid on;

% 3. Plot dell'Errore di Posizione
error_pos = xi_sim - xi_d_all;
figure('Name', 'Errore di Inseguimento', 'Color', 'w');
plot(t, error_pos(1,:), 'r', t, error_pos(2,:), 'g', t, error_pos(3,:), 'c');
grid on;
xlabel('Tempo [s]'); ylabel('Errore [m]');
legend('Errore x', 'Errore y', 'Errore z');
title('Evoluzione dell''errore nel tempo');



