%% ARS5 Project
% Title:    Quadcopter Backstepping Controller
% 
% Authors:  Stefano MUSSO PIANTELLI, 
%           Elias Charbel SALAMEH
% 
% Date:     15/01/2026

clear; 
clc; 
close all;

%% Parameters
% Physical parameters
g       = 9.81;
m       = 0.6;
Ixx     = 0.0014; Iyy = 0.0165; Izz = 0.0152;
I       = diag([Ixx Iyy Izz]);
F_ext   = [-0.0*randn(1,1);-0.1*randn(1,1);-g];

% Simulation parameters
Tend    = 40;
dt      = 1e-3;        % smaller dt for better attitude dynamics resolution
t       = 0:dt:Tend;
N       = numel(t);

%% Desired trajectory – Ascending spiral
%run('Trajectory_spiral.m');
run('Trajectory_lemniscate.m');
%run('Trajectory_cylindrical.m');
%run('Target_point.m');

%% State vector: X = [xi(3); v(3); q(4); Omega(3)]
% quaternion q = [qw; qx; qy; qz], unit norm
X       = zeros(13, N);

% Initial conditions:
xi0     = [0; 0; 0];              % initial position
%xi0     = xi_d(0);                  % instead of [0;0;0]
v0      = [0; 0; 0];                % initial velocity
q0      = [1; 0; 0; 0];             % initial attitude (no rotation)
Omega0  = [0; 0; 0];                % initial angular velocity
X(:,1)  = [xi0; v0; q0; Omega0];

%% Gains
% k1  = rand(1);
% k2  = rand(1);
% k3  = rand(1);
% k4  = rand(1);
k1  = 1;
k2  = 1;
k3  = 1;
k4  = 1;
a   = 1;
d   = 0.0;

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
tau = zeros(3,N);
F_des_norm = zeros(1,N);

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

    % variable yaw logic
    % Calculate desired heading based on planar velocity
    vx_des = xi_d_dot_val(1);
    vy_des = xi_d_dot_val(2);
    horizontal_speed = sqrt(vx_des^2 + vy_des^2);

    if horizontal_speed > 0.1 % Threshold (0.1 m/s)
        % if moving: look where we are going
        psi_d(k) = atan2(vy_des, vx_des);
    else
        % if hovering or vertical: ignore new 
        if k > 1
            psi_d(k) = psi_d(k-1);
        else
            psi_d(k) = 0; % init
        end
    end
    % --- 

    e1(:,k) = -(a*(xi_k - xi_d_val)-d); % e1
    xi_v_dot(:,k) = (k1*e1(:,k))/a + xi_d_dot_val;
    e2(:,k) = xi_dot_k - xi_v_dot(:,k); % e2
    e1_dot(:,k) = -a*e2(:,k)-k1*e1(:,k); % e1dot
    % F_des: control input
    vec1 = m * ( (a * e1(:,k)) - F_ext + ...
                 ((k1 * e1_dot(:,k)) / a + xi_d_ddot_val) - ...
                 (k2 * e2(:,k)) );

    q_row = q_k.'; %q as a row 
    q_left = quatinv(q_row);
    vec_quat = [0, vec1'];

    %compute F_u
    F_des_quat = quatmultiply(quatmultiply(q_left, vec_quat), q_row);
    F_des = F_des_quat(2:4)';

    Fu(:,k) = F_des;
    F_des_norm(:,k) = norm(F_des);

    % e2_dot
    vec2 = Fu(:,k)/m;
    q_conj2 = quatconj(q_row);
    vec_quat2 = [0, vec2'];
    e2_dot_quat = quatmultiply(quatmultiply(q_row, vec_quat2), q_conj2) ...
                    + F_ext - (k1*e1_dot(:,k)/a + xi_d_ddot_val); 
    e2_dot(:,k) = e2_dot_quat(2:4); % e2dot

    % attitude extraction
    if F_des_norm(:,k) < 1e-6
        % Fallback: Zero thrust and identity orientation (or keep previous)
        Fu_des_norm = 0;
        q_d = [1, 0, 0, 0]; % Identity quaternion [w x y z]
        warning('Fu magnitude is near zero. Thrust set to 0, attitude set to identity.');
        % tau = 0;
    end
    
    z_d = F_des / F_des_norm(:,k);
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
    q_bar_e  = quatmultiply(q_d_conj, q_k');
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
    
    vec_rotated = quatmultiply(quatmultiply(q_row, vec_quat2), q_conj2);  % 1x4
    xi_ddot = vec_rotated(2:4).' + F_ext;  % 3x1
    q_dot = 0.5*quatmultiply(q_k',[0 (omega_k.')]);

    omega_dot = I\(tau-cross(omega_k,I*omega_k));

    % noise
    noise_pos   = 0.001 * randn(3,1);    % Position noise (usually 0 for physics, this is mostly sensor noise)
    noise_vel   = 0.001 * randn(3,1);    % Velocity noise (Wind gusts, drag irregularities) - in m/s
    noise_q     = 0.001 * randn(4,1);    % Attitude noise (Vibrations affecting integration)
    noise_omega = 0.001 * randn(3,1);    % Angular velocity noise (Torque disturbances) - in rad/s

    % update
    X(1:3,k+1)   = xi_k       + dt * xi_dot + noise_pos;
    X(4:6,k+1)   = xi_dot_k   + dt * xi_ddot + noise_vel;
    
    q_next       = q_k        + dt * q_dot.' + noise_q;
    X(7:10,k+1)  = q_next / norm(q_next);
    
    X(11:13,k+1) = omega_k    + dt * omega_dot(:, k) + noise_omega;
end

%% computation of the desired position xi
xi_d_all = zeros(3, N);
for k = 1:N
    xi_d_all(:,k) = xi_d(t(k));
end

%% Plots
xi_sim = X(1:3, :);

% 3D position evolution
figure('Name', '3D Position', 'Color', 'w');
plot3(xi_d_all(1,:), xi_d_all(2,:), xi_d_all(3,:), 'g--', 'LineWidth', 1.5); hold on;
plot3(xi_sim(1,:), xi_sim(2,:), xi_sim(3,:), 'r', 'LineWidth', 1.5);
grid on; axis equal;
xlabel('X [m]'); ylabel('Y [m]'); zlabel('Z [m]');
legend('Desired (Spiral)', 'Simulated (Backstepping)');
title('3D Trajectory Tracking');
view(45, 30);

% Temporal analysis for each axis (x, y, z)
figure('Name', 'Axis Analysis x, y, z', 'Color', 'w');
subplot(3,1,1);
plot(t, xi_d_all(1,:), 'g--', t, xi_sim(1,:), 'r');
ylabel('x [m]'); grid on;
legend('Desired', 'Simulated');
title('Position Tracking on individual axes');

subplot(3,1,2);
plot(t, xi_d_all(2,:), 'g--', t, xi_sim(2,:), 'r');
ylabel('y [m]'); grid on;

subplot(3,1,3);
plot(t, xi_d_all(3,:), 'g--', t, xi_sim(3,:), 'r');
ylabel('z [m]'); xlabel('Time [s]'); grid on;

% Position Error Plot
error_pos = xi_sim - xi_d_all;
figure('Name', 'Tracking Error', 'Color', 'w');
plot(t, error_pos(1,:), 'r', t, error_pos(2,:), 'g', t, error_pos(3,:), 'b');
grid on;
xlabel('Time [s]'); ylabel('Error [m]');
legend('Error x', 'Error y', 'Error z');
title('Error evolution over time');

% control input applied
figure('Name', 'Control Inputs', 'Color', 'w');
subplot(4,1,1);
plot(t, F_des_norm(1,:), 'b');
ylabel('Thrust'); grid on;
title('Thrust w.r.t. time');

subplot(4,1,2);
plot(t, Fu(1,:), 'b');
ylabel('roll torque'); grid on;

subplot(4,1,3);
plot(t, Fu(2,:), 'b');
ylabel('pitch torque'); grid on;

subplot(4,1,4);
plot(t, Fu(3,:), 'b');
ylabel('yaw torque'); grid on;

% %% 3D Animation
% fprintf('Starting 3D Animation...\n');
% 
% % Animation Parameters
% anim_speed = 100; % Skip frames to make animation faster (Plot every 100th step)
% scale_arrow = 0.5; % Length of the direction arrow
% 
% % Create Figure
% figure('Name', '3D Drone Animation', 'Color', 'w');
% axis equal;
% grid on;
% hold on;
% view(45, 30);
% xlabel('X [m]'); ylabel('Y [m]'); zlabel('Z [m]');
% title('Quadcopter Trajectory Animation');
% 
% % limits (adjust based on your trajectory to keep view stable)
% xlim([min(xi_d_all(1,:))-2, max(xi_d_all(1,:))+2]);
% ylim([min(xi_d_all(2,:))-2, max(xi_d_all(2,:))+2]);
% zlim([min(xi_d_all(3,:))-2, max(xi_d_all(3,:))+2]);
% 
% % Draw the desired trajectory (Static background)
% plot3(xi_d_all(1,:), xi_d_all(2,:), xi_d_all(3,:), 'g--', 'LineWidth', 1, 'DisplayName', 'Desired Path');
% 
% % Initialize graphics objects (Drone, Trail, Arrow)
% % Drone body (represented as a blue dot)
% h_drone = plot3(0, 0, 0, 'bo', 'MarkerFaceColor', 'b', 'MarkerSize', 8, 'DisplayName', 'Drone');
% 
% % Trail (path taken so far)
% h_trail = plot3(0, 0, 0, 'r-', 'LineWidth', 1, 'DisplayName', 'Flown Path');
% 
% % Orientation Arrow (Red line indicating Body-X / Heading)
% h_arrow = line([0 0], [0 0], [0 0], 'Color', 'y', 'LineWidth', 2, 'DisplayName', 'Heading');
% 
% legend('show');
% 
% % Animation Loop
% for k = 1:anim_speed:N
%     % Extract current state
%     pos = X(1:3, k);       % Position [x; y; z]
%     quat = X(7:10, k);     % Quaternion [w; x; y; z]
% 
%     % Update Drone Position
%     set(h_drone, 'XData', pos(1), 'YData', pos(2), 'ZData', pos(3));
% 
%     % Update Trail
%     % (For performance, update trail in chunks or just the current data)
%     set(h_trail, 'XData', X(1, 1:k), 'YData', X(2, 1:k), 'ZData', X(3, 1:k));
% 
%     % Calculate Heading Arrow
%     % Convert quaternion to Rotation Matrix
%     % Note: quat2rotm expects [w x y z] 
%     R_current = quat2rotm(quat'); 
% 
%     % The Body-X axis in World Frame is the first column of R
%     heading_vec = R_current(:, 1); 
% 
%     % Calculate arrow end point
%     arrow_end = pos + (heading_vec * scale_arrow);
% 
%     % Update Arrow
%     set(h_arrow, 'XData', [pos(1), arrow_end(1)], ...
%         'YData', [pos(2), arrow_end(2)], ...
%         'ZData', [pos(3), arrow_end(3)]);
% 
%     % Update Title with time
%     title(sprintf('Simulation Time: %.2f s', t(k)));
% 
%     % Force draw
%     drawnow limitrate; 
%     pause(0.01);
% end


%% 3D Animation + Video Export
fprintf('Starting 3D Animation...\n');

% Animation Parameters
anim_speed  = 100;      % Plot every anim_speed steps
scale_arrow = 0.5;      % Length of heading arrow

% ===================== FIGURE =====================
figure('Name', '3D Drone Animation', 'Color', 'w');
axis equal;
grid on;
hold on;
view(45, 30);
xlabel('X [m]'); ylabel('Y [m]'); zlabel('Z [m]');
title('Quadcopter Trajectory Animation');

% Limits (keep camera stable)
xlim([min(xi_d_all(1,:))-2, max(xi_d_all(1,:))+2]);
ylim([min(xi_d_all(2,:))-2, max(xi_d_all(2,:))+2]);
zlim([min(xi_d_all(3,:))-2, max(xi_d_all(3,:))+2]);

% Desired trajectory (background)
plot3(xi_d_all(1,:), xi_d_all(2,:), xi_d_all(3,:), ...
      'g--', 'LineWidth', 1, 'DisplayName', 'Desired Path');

% ===================== GRAPHICS OBJECTS =====================
% Drone body
h_drone = plot3(0, 0, 0, 'bo', ...
    'MarkerFaceColor', 'b', 'MarkerSize', 8, 'DisplayName', 'Drone');

% Trail
h_trail = plot3(0, 0, 0, 'r-', ...
    'LineWidth', 1, 'DisplayName', 'Flown Path');

% Heading arrow (Body-X axis)
h_arrow = line([0 0], [0 0], [0 0], ...
    'Color', 'y', 'LineWidth', 2, 'DisplayName', 'Heading');

legend('show');

% ===================== VIDEO WRITER =====================
video_name = 'Quadcopter_Animation.mp4';
v = VideoWriter(video_name, 'MPEG-4');
v.FrameRate = 30;     % frames per second
v.Quality   = 100;    % max quality
open(v);

% ===================== ANIMATION LOOP =====================
for k = 1:anim_speed:N

    % Current state
    pos  = X(1:3, k);
    quat = X(7:10, k);   % [w x y z]

    % Drone position
    set(h_drone, 'XData', pos(1), 'YData', pos(2), 'ZData', pos(3));

    % Trail
    set(h_trail, 'XData', X(1,1:k), ...
                 'YData', X(2,1:k), ...
                 'ZData', X(3,1:k));

    % Heading arrow
    R_current  = quat2rotm(quat');      % Rotation matrix
    heading_vec = R_current(:,1);       % Body-X in world
    arrow_end  = pos + scale_arrow * heading_vec;

    set(h_arrow, 'XData', [pos(1), arrow_end(1)], ...
                 'YData', [pos(2), arrow_end(2)], ...
                 'ZData', [pos(3), arrow_end(3)]);

    % Title with time
    title(sprintf('Simulation Time: %.2f s', t(k)));

    % Render & save frame
    drawnow limitrate;
    frame = getframe(gcf);
    writeVideo(v, frame);

end

% ===================== CLOSE VIDEO =====================
close(v);
fprintf('Video successfully saved: %s\n', video_name);
