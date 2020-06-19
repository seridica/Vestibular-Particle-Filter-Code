%%%
% File: MacNeilageParticleModel.m
% Author: Calvin Kuo
% Date: 11-15-2018
% Notes: Runs the MacNeilage Particle model equivalent

%% Model parameters
% Parameters
%   - Canal Dynamics
%       : tc = Canal time constant

%dt = 0.0005;
dt = 1/8000;
Fs = 1/dt;

N = 100;

% Parameters from Borah Paper
params.Q = ones(1,N) * 200000000;
params.R = ones(1,N) * 25;
params.tc = ones(1,N) * 1/4;
params.beta = ones(1,N) * 30;
params.N = N;

%% Inputs
% Motion profile (angular acceleration)

% Motion Profile
mopo = 2;
[t, angAcc] = MotionProfile( mopo, dt, 60.5 );
%angAcc = cumtrapz( angAcc(2,:)' ) * dt;
angAcc = angAcc(2,:)';

%% Simulating input with noise and filtered according to bandwidth

%% Laurens particle filter
AllStates = zeros(5, length(t));
x_curr = zeros(4,N);
x_curr(4,:) = normrnd( 0, sqrt(params.R) );
for i=1:length( t )
    if ( i*dt < 1.0 )
        st = 0;
    else
        st = 1;
    end
    x_new = MacNeilageParticleDeriv( x_curr, angAcc(i), params, dt, st );
    
    AllStates(:,i) = mean( x_new' )';
    x_curr = x_new(1:4,:);
end