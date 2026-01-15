%% ARS5 Project
% Title:    Multi-Agent Platooning (Daisy-Chain)
% Authors:  Stefano MUSSO PIANTELLI, Elias Charbel SALAMEH
% Date:     15/01/2026

clear; clc; close all;

%% Parameters
num_agents = 3;

g = 9.81; 
m = 0.6;
I = diag([0.0014, 0.0165, 0.0152]);

Tend = 40;
dt   = 1e-3;
t    = 0:dt:Tend;
N    = numel(t);

%% Desired trajectory – Ascending spiral
% run('Trajectory_spiral.m');
run('Trajectory_lemniscate.m');
% run('Trajectory_cylindrical.m');
% run('Target_point.m');

%% State Initialization
X = zeros(13, N, num_agents); 

% Robot 1: Leader
xi0_1 = [10; -5; 7];
% Robot 2: 2m behind Leader
xi0_2 = xi0_1 + [-2; 0; 0];
% Robot 3: 2m behind Robot 2
xi0_3 = xi0_2 + [-2; 0; 0];

xi0_all = [xi0_1, xi0_2, xi0_3];

for i = 1:num_agents
    X(:, 1, i) = [xi0_all(:,i); zeros(3,1); [1;0;0;0]; zeros(3,1)];
end

%% Platooning Topology (1 -> 2 -> 3)
A = zeros(num_agents, num_agents);
A(2,1) = 1; % Robot 2 listens to Robot 1
A(3,2) = 1; % Robot 3 listens to Robot 2 (NOT Robot 1)

% Desired Relative Offsets (xi_i - xi_j)_desired
% "Stay 2m behind the one you are following"
D_offset = zeros(3, num_agents, num_agents);
D_offset(:, 2, 1) = [-2; -1.5; -0.3]; % R2 relative to R1
D_offset(:, 3, 2) = [-2; -0.4; -0.3]; % R3 relative to R2

%% Gains
k1 = 1; 
k2 = 1; 
k3 = 1; 
k4 = 1;
a  = 1;

%% Pre-allocation
e1     = zeros(3, N, num_agents);
e1_dot = zeros(3, N, num_agents);
e2     = zeros(3, N, num_agents);
Fu     = zeros(3, N, num_agents);
psi_d  = zeros(1, N, num_agents);

%% Main Loop
for k = 1:N-1
    tk = t(k);
    
    % Global Trajectory (Leader reference)
    xi_d_val      = xi_d(tk);
    xi_d_dot_val  = xi_d_dot(tk);
    xi_d_ddot_val = xi_d_ddot(tk);
    
    for i = 1:num_agents
        
        % State Extraction
        xi_i     = X(1:3, k, i);
        v_i      = X(4:6, k, i);
        q_i      = X(7:10, k, i);
        omega_i  = X(11:13, k, i);
        
        % Consensus Error
        sum_pos_error = [0;0;0];
        sum_vel_ref   = [0;0;0];
        sum_acc_ref   = [0;0;0];
        degree_connectivity = 0;
        
        % R1
        if i == 1
            sum_pos_error = (xi_i - xi_d_val);
            sum_vel_ref   = xi_d_dot_val;
            sum_acc_ref   = xi_d_ddot_val;
            degree_connectivity = 1;
        
        % R2, 3
        else
            % Iterate to find who I am following
            for j = 1:num_agents
                if A(i,j) == 1
                    xi_j = X(1:3, k, j);
                    v_j  = X(4:6, k, j);
                    
                    % Error: xi_i-xi_k-d_ij
                    sum_pos_error = sum_pos_error + (xi_i - xi_j - D_offset(:, i, j));
                    
                    % Feedforward Velocity: Copy his velocity
                    sum_vel_ref   = sum_vel_ref + v_j; 
                    
                    degree_connectivity = degree_connectivity + 1;
                end
            end
        end
        
        e1(:, k, i) = - (a * sum_pos_error);
        
        xi_v_dot = (sum_vel_ref + k1 * e1(:, k, i)) / max(1, degree_connectivity); 
        
        e2(:, k, i) = v_i - xi_v_dot;
        e1_dot(:, k, i) = -a * e2(:, k, i) - k1 * e1(:, k, i);
        
        F_ext = [-0.1*randn;-0.1*randn;-g]; 
        
        vec1 = m * ( (a * e1(:,k,i)) - F_ext + ...
                     ((k1 * e1_dot(:,k,i))/a + sum_acc_ref) - ...
                     (k2 * e2(:,k,i)) );
                 
        q_row  = q_i.';
        q_left = quatinv(q_row);
        vec_quat = [0, vec1'];
        F_des_quat = quatmultiply(quatmultiply(q_left, vec_quat), q_row);
        F_des = F_des_quat(2:4)';
        
        Fu(:, k, i) = F_des;
        F_des_norm = norm(F_des);
        
        % Platooning Yaw
        % Leader: Look at Velocity
        if i == 1
            vel_ref = xi_d_dot_val;
            if norm(vel_ref(1:2)) > 0.1
                psi_d(1, k, i) = atan2(vel_ref(2), vel_ref(1));
            else
                if k>1, psi_d(1, k, i) = psi_d(1, k-1, i); else, psi_d(1, k, i) = 0; end
            end
        else
            % Follower: Look at Predecessor
            % Find 'j' where A(i,j)=1
            [~, predecessor_idx] = max(A(i,:)); 
            
            % Vector pointing to predecessor
            look_vec = X(1:3, k, predecessor_idx) - xi_i;
            
            if norm(look_vec(1:2)) > 0.1
                 psi_d(1, k, i) = atan2(look_vec(2), look_vec(1));
            else
                 if k>1, psi_d(1, k, i) = psi_d(1, k-1, i); else, psi_d(1, k, i) = 0; end
            end
        end
        
        % Attitude Calculation
        z_d = F_des / max(F_des_norm, 1e-6);
        l_vec = [-sin(psi_d(1, k, i)); cos(psi_d(1, k, i)); 0];
        x_d = cross(l_vec, z_d); x_d = x_d/norm(x_d);
        y_d = cross(z_d, x_d);
        Rd = [x_d, y_d, z_d];
        q_d = rotm2quat(Rd);
        
        q_i = q_i / norm(q_i); q_d = q_d / norm(q_d);
        q_d_conj = quatconj(q_d);
        q_bar_e  = quatmultiply(q_d_conj, q_i');
        q_e = [1 - abs(q_bar_e(1)); q_bar_e(2:4).'];
        
        M = quatM(q_bar_e);
        omega_v = -2*k3*M.'*q_e;
        alpha_2 = omega_i - omega_v;
        q_e_dot = 0.5 * M * (alpha_2 + omega_v);
        M_dot = quatM(q_e_dot);
        omega_v_dot = -2*k3*(M_dot.' * q_e + M.' * q_e_dot);
        
        tau = cross(omega_i, I*omega_i) + I*omega_v_dot - k4*alpha_2 - 0.5*M.'*q_e;
        
        % Integration
        vec2 = F_des/m;
        q_conj2 = quatconj(q_row);
        vec_quat2 = [0, vec2'];
        
        vec_rotated = quatmultiply(quatmultiply(q_row, vec_quat2), q_conj2);
        xi_ddot = vec_rotated(2:4).' + F_ext;
        
        q_dot = 0.5*quatmultiply(q_i', [0, omega_i']);
        omega_dot = I \ (tau - cross(omega_i, I*omega_i));
        
        noise_pos = 0.001 * randn(3,1); noise_vel = 0.001 * randn(3,1);
        
        X(1:3, k+1, i)   = xi_i     + dt * v_i + noise_pos;
        X(4:6, k+1, i)   = v_i      + dt * xi_ddot + noise_vel;
        q_next           = q_i      + dt * q_dot.';
        X(7:10, k+1, i)  = q_next / norm(q_next);
        X(11:13, k+1, i) = omega_i  + dt * omega_dot;
    end
end

%% Visualization (Platooning)
xi_d_all = zeros(3, N);
for k = 1:N, xi_d_all(:,k) = xi_d(t(k)); end

colors = {'r', 'b', 'g'};
labels = {'Leader', 'Follower 1', 'Follower 2'};

figure('Name', 'Platoon Trajectory', 'Color', 'w');
plot3(xi_d_all(1,:), xi_d_all(2,:), xi_d_all(3,:), 'g--', 'LineWidth', 1.5, 'DisplayName', 'Reference'); 
hold on; grid on; axis equal; view(45, 30);
xlabel('X'); ylabel('Y'); zlabel('Z');

for i = 1:num_agents
    pos = squeeze(X(1:3, :, i));
    plot3(pos(1,:), pos(2,:), pos(3,:), 'Color', colors{i}, 'LineWidth', 1.5, 'DisplayName', labels{i});
    plot3(pos(1,end), pos(2,end), pos(3,end), 's', 'Color', colors{i}, 'MarkerFaceColor', colors{i});
end
legend('show');
title('Platooning (Daisy-Chain) Control');

%% 3D Animation + Video Export (MULTI-ROBOT – IMPROVED)
fprintf('Starting 3D Animation...\n');

anim_speed  = 100;
scale_arrow = 0.6;
margin      = 2;   % extra space around trajectories

%% ===== PRECOMPUTE GLOBAL AXIS LIMITS (ALL ROBOTS) =====
all_pos = reshape(X(1:3,:,:), 3, []);
xmin = min(all_pos(1,:)) - margin;
xmax = max(all_pos(1,:)) + margin;
ymin = min(all_pos(2,:)) - margin;
ymax = max(all_pos(2,:)) + margin;
zmin = min(all_pos(3,:)) - margin;
zmax = max(all_pos(3,:)) + margin;

%% ===== FIGURE =====
figure('Name', '3D Platoon Animation', 'Color', 'w');
axis equal; grid on; hold on;
view(45,30);
xlabel('X [m]'); ylabel('Y [m]'); zlabel('Z [m]');
title('Multi-Agent Platooning Animation');

xlim([xmin xmax]);
ylim([ymin ymax]);
zlim([zmin zmax]);

%% ===== REFERENCE TRAJECTORY =====
h_ref = plot3(xi_d_all(1,:), xi_d_all(2,:), xi_d_all(3,:), ...
              'k--', 'LineWidth', 1.5);

%% ===== GRAPHICS OBJECTS =====
colors  = {'r','b','g'};
labels  = {'Leader','Follower 1','Follower 2'};

h_drone = gobjects(num_agents,1);
h_trail = gobjects(num_agents,1);
h_arrow = gobjects(num_agents,1);

for i = 1:num_agents
    h_drone(i) = plot3(0,0,0,'o', ...
        'Color', colors{i}, ...
        'MarkerFaceColor', colors{i}, ...
        'MarkerSize', 8);

    h_trail(i) = plot3(0,0,0,'-', ...
        'Color', colors{i}, ...
        'LineWidth', 1.2);

    h_arrow(i) = line([0 0],[0 0],[0 0], ...
        'Color', colors{i}, ...
        'LineWidth', 2);
end

%% ===== CORRECT LEGEND (USING HANDLES) =====
legend([h_ref; h_drone], ...
       [{'Reference'}, labels], ...
       'Location','best');

%% ===== VIDEO WRITER =====
video_name = 'Platooning_Animation.mp4';
v = VideoWriter(video_name,'MPEG-4');
v.FrameRate = 30;
v.Quality   = 100;
open(v);

%% ===== ANIMATION LOOP =====
for k = 1:anim_speed:N

    for i = 1:num_agents
        pos  = X(1:3, k, i);
        quat = X(7:10, k, i);

        % Drone body
        set(h_drone(i), ...
            'XData', pos(1), ...
            'YData', pos(2), ...
            'ZData', pos(3));

        % Trail
        set(h_trail(i), ...
            'XData', X(1,1:k,i), ...
            'YData', X(2,1:k,i), ...
            'ZData', X(3,1:k,i));

        % Heading arrow (body X-axis)
        R = quat2rotm(quat');
        heading = R(:,1);
        arrow_end = pos + scale_arrow * heading;

        set(h_arrow(i), ...
            'XData', [pos(1) arrow_end(1)], ...
            'YData', [pos(2) arrow_end(2)], ...
            'ZData', [pos(3) arrow_end(3)]);
    end

    title(sprintf('Simulation Time: %.2f s', t(k)));
    drawnow limitrate;

    frame = getframe(gcf);
    writeVideo(v, frame);
end

%% ===== CLOSE VIDEO =====
close(v);
fprintf('Video successfully saved: %s\n', video_name);

