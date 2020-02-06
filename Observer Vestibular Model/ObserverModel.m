%%%
% File: ObserverModel.m
% Author: Calvin Kuo
% Date: 8-23-2018
% Notes: Runs the observer model

%% Model parameters
% Parameters
%   - Canal Dynamics
%       : tc = Canal time constant
%       : K = Observer gain

params.tc = [4;4;4];
params.K = 3;

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

%% Simulating with Observer Model - Deriv
x_prev = zeros(9,1);%[randn(12,1);zeros(3,1)];
full_states = zeros(9,length(t));

for i=2:length( t )
    alpha = angAcc(:,i);
    x_new = ObserverDeriv( x_prev, alpha, params, dt );
    
    full_states(:,i) = x_new;
    x_prev = x_new;
end

figure(1); clf; hold on;
plot( t, full_states(2,:), 'k-' );
plot( t, full_states(5,:), 'r-' );
plot( t, full_states(8,:), 'g-' );

%% Simulating with Observer Model - TRANSFER FUNCTION
figure(2); clf;

canal_tf_x = tf([1], [1,1/params.tc(1)]);
canal_tf_y = tf([1], [1,1/params.tc(2)]);
canal_tf_z = tf([1], [1,1/params.tc(3)]);

canal_model_x = tf([1,0], [1,1/params.tc(1)]);
canal_model_y = tf([1,0], [1,1/params.tc(2)]);
canal_model_z = tf([1,0], [1,1/params.tc(3)]);

[y,t,x] = lsim([canal_tf_y; canal_model_y*(params.K*canal_tf_y) / (1+params.K*canal_model_y); (params.K*canal_tf_y) / (1+params.K*canal_model_y)], angAcc(2,:), t);
clf; hold on;
plot( t, y(:,1), 'k-' );
plot( t, y(:,2), 'r-' );
plot( t, y(:,3), 'g-' )
xlim([-50,150])