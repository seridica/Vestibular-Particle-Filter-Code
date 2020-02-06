%%%
% File: RaphanAndCohenModel.m
% Author: Calvin Kuo
% Date: 8-23-2018
% Notes: Runs the most basic Raphan and Cohen velocity storage model.

%% Model parameters
% Parameters
%   - Canal Dynamics
%       : tc = Canal time constant
%       : tvs = Velocity leakage time constant

params.tc = [4;4;4];
params.tvs = [10;10;10];
params.g = 0.25;

dt = 0.1;
Fs = 1/dt;

%% Inputs
% Motion profile (angular acceleration)

% Motion Profile
mopo = 2;
[t, angAcc] = MotionProfile( mopo, dt );

% Angular position (for now)
params.Tcan = [1 0 0;...
               0 1 0;...
               0 0 1];

% params.Tcan = [1 0 0; ...
%                0 cos(pi/4) -sin(pi/4); ...
%                0 sin(pi/4) cos(pi/4)];
% 
% params.Tcan = [1 0 0; ...
%                0 cos(pi/4) -sin(pi/4); ...
%                0 sin(pi/4) cos(pi/4)] * ...
%               [cos(pi/4) -sin(pi/4) 0; ...
%                sin(pi/4) cos(pi/4) 0; ...
%                0 0 1];
% params.Tcan = [cos(pi/2) -sin(pi/2) 0; ...
%                sin(pi/2) cos(pi/2) 0; ...
%                0 0 1];

%% Simulating with Raphan and Cohen - DERIV
x_prev = zeros(6,1);%[randn(12,1);zeros(3,1)];
full_states = zeros(6,length(t));

for i=2:length( t )
    alpha = angAcc(:,i);
    x_new = RaphanAndCohenDeriv( x_prev, alpha, params, dt );
    
    full_states(:,i) = x_new;
    x_prev = x_new;
end

figure(1); clf; hold on;
plot( t, full_states(2,:), 'k-' );
plot( t, full_states(5,:), 'r-' );
plot( t, full_states(2,:) + full_states(5,:), 'g-' );

%% Simulating with Raphan and Cohen - TRANSFER FUNCTION
figure(2); clf;

canal_tf_x = tf([1], [1,1/params.tc(1)]);
canal_tf_y = tf([1], [1,1/params.tc(2)]);
canal_tf_z = tf([1], [1,1/params.tc(3)]);

velstore_tf_x = tf([1], [1,1/params.tvs(1)]);
velstore_tf_y = tf([1], [1,1/params.tvs(2)]);
velstore_tf_z = tf([1], [1,1/params.tvs(3)]);

[y,t,x] = lsim([canal_tf_y; params.g*canal_tf_y * velstore_tf_y], angAcc(2,:), t);
clf; hold on;
plot( t, y(:,1), 'k-' );
plot( t, y(:,2), 'r-' );
plot( t, y(:,1) + y(:,2), 'g-' )
xlim([-50,150])