%%%
% File: MacNeilageKalmanBucyModel.m
% Author: Calvin Kuo
% Date: 11-12-2018
% Notes: Runs the Laurens Particle model

%% Model parameters
% Parameters
%   - Canal Dynamics
%       : tc = Canal time constant

dt = 0.001;
Fs = 1/dt;

n_states = 13;

% Parameters from Borah Paper
params.Q = 200000000; %1000000000;
params.R = 25;
params.tc = 1/4;
params.beta = 30;

%% Inputs
% Motion profile (angular acceleration)

% Motion Profile
mopo = 2;
[t, angAcc] = MotionProfile( mopo, dt, 60.5 );
%angAcc = cumtrapz( angAcc(2,:)' ) * dt;
angAcc = angAcc(2,:)';

%% Simulating input with noise and filtered according to bandwidth

%% Laurens particle filter
AllStates = zeros(n_states, length(t));
x_curr = [zeros(4,1); reshape( eye(3), 9, 1 ) * 0.001];
for i=1:length( t )
    x_new = MacNeilageKalmanBucyDeriv( x_curr, angAcc(i), params, dt );
    
    AllStates(:,i) = x_new;
    x_curr = x_new;
end