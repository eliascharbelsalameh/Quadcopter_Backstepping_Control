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
Tend    = 40;
dt      = 1e-3;        % smaller dt for better attitude dynamics resolution
t       = 0:dt:Tend;
N       = numel(t);

%% Desired trajectory – Ascending spiral
R      = 1.0;          % [m] spiral radius
omega_traj  = 0.5;     % [rad/s] angular speed
z0     = 0.0;          % [m] initial altitude
zf     = 1.5;          % [m] final altitude at t = Tend
vz     = (zf - z0)/Tend;  % [m/s] vertical speed

% Position
xi_d = @(tt) [ R*cos(omega_traj*tt);         % x
    R*sin(omega_traj*tt);         % y
    z0 + vz*tt ];            % z

% Velocity
xi_d_dot = @(tt) [ -R*omega_traj*sin(omega_traj*tt);   % ẋ
    R*omega_traj*cos(omega_traj*tt);   % ẏ
    vz ];                    % ż (constant)

% Acceleration
xi_d_ddot = @(tt) [ -R*omega_traj^2*cos(omega_traj*tt);  % ẍ
    -R*omega_traj^2*sin(omega_traj*tt);  % ÿ
    0 ];                      % z̈

%% State vector: X = [xi(3); v(3); q(4); Omega(3)]
% quaternion q = [qw; qx; qy; qz], unit norm
X       = zeros(13, N);

% Initial conditions:
% xi0     = [0; 0; 0];                % initial position
xi0     = xi_d(0);                  % instead of [0;0;0]
v0      = [0; 0; 0];                % initial velocity
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
e2 = zeros(3, N);
qe = zeros(4, N);
alpha_2 = zeros(3, N);

xi_v_dot = zeros(3,N);
xi_v_ddot = zeros(3,N);
omega_v = zeros(3,N);
omega_v_dot = zeros(3,N);
qe_dot = zeros(4,N);

%% sampled loop
for k = 1:N-1
    tk = t(k);

    xi_k = x(1:3,k);
    xi_dot_k = x(4:6,k);
    q_k = x(7:10,k);
    omega_k = x(11:13,k);

    % target
    xi_d      = xi_d(tk);
    xi_d_dot  = xi_d_dot(tk);
    xi_d_ddot = xi_d_ddot(tk);  

    e1(k) = -(a*(xi_k - xi_d_dot)-d);
    xi_v_dot(k) = (k1*e1(k))/a + xi_d_dot; 
    e2(k) = xi_k - xi_v_dot;
    e1_dot = -a*e2(k)-k1*e1(k);
    % F_des: control input
    vec1 = m * ( (a * e1) - F_ext(tk) + ...
                 ((k1 * e1_dot) / a + xi_target_ddot) - ...
                 (k2 * e2) );
    q_row = q'; 
    q_left = quatinv(q_row);
    q_conj = quatconj(q_row);
    q_right = quatinv(q_conj);
    vec_quat = [0, vec1']; 
    F_des_quat = quatmultiply(quatmultiply(q_left, vec_quat), q_right);
    F_des = F_des_quat(2:4)';
    % Fu
    Fu = norm(F_des);
    % e2_dot
    vec2 = F_u/m;
    q_conj2 = quatconj(q_row);
    vec_quat2 = [0, vec2'];
    e2_dot_quat = quatmultiply(quatmultiply(q_row, vec_quat2), q_conj2) + F_ext(tk) - (k1*e1_dot/a + xi_d_ddot);
    e2_dot = e2_dot_quat(2:4);
    % attitude extraction
    if Fu < 1e-6
        % Fallback: Zero thrust and identity orientation (or keep previous)
        Fu = 0;
        q_d = [1, 0, 0, 0]; % Identity quaternion [w x y z]
        warning('U_pos magnitude is near zero. Thrust set to 0, attitude set to identity.');
        return;
    end
    z_d = U_pos / Fu;
    l_vec = [-sin(psi_d); cos(psi_d); 0];

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
    q_bar_e  = quatmultiply(q_d_conj, q_k);
    % q_e from q_bar_e
    q_bar_e0 = q_bar_e(1); q_bar_ev = q_bar_e(2:4).';
    q_e = [1 - abs(q_bar_e0); q_bar_ev];
    % M : quaternion error dynamics computation
    M = quatM(q_bar_e);
    omega_v = -2*k3*M.'*q_e;
    alpha_2 = omega_k - omega_v;
    q_e_dot = 0.5 * M * (alpha_2 + omega_v);
    M_dot = quatM(q_e_dot);
    omega_v_dot = -2*k3*(M_dot.' * q_e + M.' * q_e_dot);
    % tau
    tau = cross(omega_k,eye(3)*omega_k)+eye(3)*omega_v_dot-k4*alpha_2-0.5*M*q_e.';
    alpha_2_dot = I\(-cross(omega_k,eye(3)*omega_k)+tau)-omega_v_dot;

    % initial system
    xi_dot = xi_dot_k;
    
    q_k = quatmultiply(quatinv(q_conj),q_e);
    q_k = q_k / norm(q_k);
    xi_ddot = quatmultiply(quatmultiply(q_row, vec_quat2), q_conj2) + F_ext(tk);
    
    q_dot = 0.5*quatmultiply(q_k,[0 (omega_k.')]);

    omega_dot = I\(tau-cross(omega_k,I*omega_k));

    % update
    x(1:3,k+1) = xi_k       + dt * xi_dot;
    x(4:6,k+1) = xi_dot_k   + dt * xi_ddot;
    x(7:10,k+1) = q_k       + dt * q_dot;
    x(11:13,k+1) = omega_k  + dt * omega_dot;
end