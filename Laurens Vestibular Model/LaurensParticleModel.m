%%%
% File: LaurensParticleModel.m
% Author: Calvin Kuo
% Date: 11-3-2018
% Notes: Runs the Laurens Particle model

%% Model parameters
% Parameters
%   - Canal Dynamics
%       : tc = Canal time constant

% From Laurens section 2.6
dt = 0.1; %1/8000; %0.1;
Fs = 1/dt;

% Parameters from Karmali Paper
%N = 3000;
N = 100;
params.N = N;
%params.tc = ones(N,1) * 4; %randn(N,1)*0.5 + 4; %ones(N,1) * 4;
params.tc = 4;
params.eta = ones(N,1) * 5;
params.omv = ones(N,1) * 10;

%% Inputs
% Motion profile (angular acceleration)

% Motion Profile
mopo = 2;
[t, angAcc] = MotionProfile( mopo, dt, 60.5 );
angAcc = angAcc(2,:)';

%% Simulating true canal dynamics
%Vfull = zeros(N, length(t));
Vfull = zeros(1, length(t));

for i=2:length(t)
    %Vfull(:,i) = ( ( -1 ./ params.tc ) .* Vfull(:,i-1) - angAcc(i) ) .* dt + Vfull(:,i-1);
    Vfull(i) = ( ( -1 ./ params.tc ) .* Vfull(i-1) + angAcc(i) ) .* dt + Vfull(i-1);
end

%% Laurens particle filter
%HeadAngVelParticles = zeros(N, length(t));
HeadAngVelParticles = zeros(1, length(t));
x_curr = zeros(N,2);
for i=1:length( t )
    V = Vfull(:,i);
    x_new = LaurensParticleDeriv( x_curr, V, params, dt );
    
    %HeadAngVelParticles(:,i) = x_new(:,2);
    HeadAngVelParticles(i) = mean( x_new(:,2) );
    x_curr = x_new;
end