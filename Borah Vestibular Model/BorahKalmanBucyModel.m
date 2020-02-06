%%%
% File: BorahKalmanBucyModel.m
% Author: Calvin Kuo
% Date: 11-12-2018
% Notes: Runs the Laurens Particle model

%% Model parameters
% Parameters
%   - Canal Dynamics
%       : tc = Canal time constant

dt = 0.001;
Fs = 1/dt;

n_states = 40;

% Parameters from Borah Paper
params.Q = 3 * 10^5;

%% Inputs
% Motion profile (angular acceleration)

% Motion Profile
mopo = 2;
[t, angAcc] = MotionProfile( mopo, dt );
angAcc = cumtrapz( angAcc(2,:)' );

%% Simulating input with noise and filtered according to bandwidth

%% Laurens particle filter
AllStates = zeros(n_states, length(t));
x_curr = [zeros(8,1); reshape( eye(4), 16, 1 ) * 0.001; reshape( eye(4), 16, 1 ) * 0.001];
for i=1:length( t )
    x_new = BorahKalmanBucyDeriv( x_curr, angAcc(i), params, dt );
    
    AllStates(:,i) = x_new;
    x_curr = x_new;
end