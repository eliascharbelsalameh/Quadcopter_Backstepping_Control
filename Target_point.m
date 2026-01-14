%% point_target.m
% Traiettoria costante verso un punto fisso

% Punto fisso desiderato
target_point = [10.0; 10.0; 10.5];  % x, y, z

% Funzioni con le stesse variabili degli altri script
xi_d = @(tt) target_point;       % posizione sempre uguale
xi_d_dot = @(tt) [0;0;0];        % velocità nulla
xi_d_ddot = @(tt) [0;0;0];       % accelerazione nulla
