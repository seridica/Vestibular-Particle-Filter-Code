%%%
% File: KarmaliKalmanBucyModel.m
% Author: Calvin Kuo
% Date: 11-13-2019
% Notes: Runs the Karmali Kalman Bucy model

%% Model parameters
% Parameters
%   - Canal Dynamics
%       : tc = Canal time constant

% From Karmali
dt = 1/8000;
%dt = 0.001; %1/8000;
Fs = 1/dt;

% Parameters from Karmali Paper
params.tc = 5.7;
params.R = 10.0;
params.Q = 41;
params.t2 = 0.005;

%% Inputs
% Motion profile (angular acceleration)

% Motion Profile
mopo = 2;
[t, angAcc] = MotionProfile( mopo, dt, 60.5 );
angAcc = cumtrapz( angAcc(2,:) ) * dt;

%% Simulating with Observer Model - Deriv
x_prev = [zeros(4,1); 0;0;0;0.001];
full_states = zeros( 12, length(t) );
    
for i=2:length( t )
    alpha = angAcc(i);
    x_new = KarmaliKalmanBucyDeriv( x_prev, alpha, params, dt );
    
    full_states(:,i) = x_new;
    x_prev = x_new(1:9);
end