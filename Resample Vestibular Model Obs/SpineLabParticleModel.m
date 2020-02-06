%%%
% File: ParticleModel.m
% Author: Calvin Kuo
% Date: 10-23-2018
% Notes: This code is a restructuring of the Laurens_NatureNeuroscience_M.m
% code originally written by PAF, CJD, NK, and JSB

%% Model parameters
% Parameters
%   - Canal Dynamics
%       : tc = Canal time constant
%       : tn = Long term adaptation time constant
%       : kv = Gain on canal afferent
%       : tvs = Velocity leakage time constant
%   - Gravity Estimation
%       : kf = Gain on gravity correction (set to 0 to turn off)
%       : ts = Gravity estimate time constant
%   - Visual Feedback
%       : go = Direct visual pathway gain (set to 0 to turn off)d
%       : ko = Indirect visual pathway gain (set to 0 to turn off)
%   - Velocity Estimate
%       : gd = Gain on canal afferent

nParticles = 50;
params.tc = ones( nParticles, 1 )* 4; %normrnd( 4.5, 0.8, [nParticles, 1] );

% For dual line
params.tn_list = [ ones(nParticles/2,1)*10000; ones(nParticles/2,1)*10000 ];

dt = 0.1;
Fs = 1/dt;

% Nois Histories
sigC = eye(3)*1.7^2; %[3.6];%[0.1]*pi/180 / dt;             % Canal noise (Laurens and Angelaki)
sigQ1 = eye(3)*1.7^2;%6.4^2;%6.4^2;%6.4^2; %[14.0];
sigQ2 = eye(3)*5.0^2;%14^2;%14^2;%6.4^2; %[14.0];

params.sigAlpha = sigC;
params.sigPrior = eye(3)*10000^2;
params.sigAfferent_list = zeros(3,3,nParticles);

for i=1:75
    params.sigAfferent_list(:,:,i) = sigQ1;
end
for i=76:100
    params.sigAfferent_list(:,:,i) = sigQ1;
end
for i=101:175
    params.sigAfferent_list(:,:,i) = sigQ1;
end
for i=176:200
    params.sigAfferent_list(:,:,i) = sigQ1;
end

%% Inputs
%   - Angular position of canals (constant 3x1 or time varying 4xn)
%   - Angular velocity (constant 3x1 or time varying 4xn)
%   - Angular acceleration (constant 3x1 or time varying 4xn)
%   - Gravity (constant 3x1 or time varying 4xn)
%   - Linear acceleration (constant 3x1 or time varying 4xn)
%   - Canal Noise (constant 3x1 or time varying 4xn)
%   - Endolymph Noise Derivative (constant 3x1 or time varying 4xn)
%   - Visual Noise (constant 3x1 or time varying 4xn)
%   - Canal orientation (constant 3x3 or time varying 4x3xn)


% Motion Profile
mopo = 2;
[t, angAcc] = MotionProfile( mopo, dt );

% Gravity
u.grav = [0; -9.81; 0];

%Angular position (for now)
params.Tcan = [1 0 0;...
               0 1 0;...
               0 0 1];

% params.Tcan = [1 0 0; ...
%                0 cos(pi/4) -sin(pi/4); ...
%                0 sin(pi/4) cos(pi/4)];

% params.Tcan = [cos(pi/2) -sin(pi/2) 0; ...
%                sin(pi/2) cos(pi/2) 0; ...
%                0 0 1];

% params.Tcan = [1 0 0; ...
%                0 cos(pi/4) -sin(pi/4); ...
%                0 sin(pi/4) cos(pi/4)] * ...
%               [cos(pi/4) -sin(pi/4) 0; ...
%                sin(pi/4) cos(pi/4) 0; ...
%                0 0 1];

%% Initializing the variables
% Initializing initial conditions
C = zeros(3,1);              % Canal output signal   
D = zeros(3,1);              % Cana afferent signal
INT = zeros(3,1);            % Integral output (Endolymph)
VS = zeros(3,1);             % Velocity storage leaky integrator
VSf = zeros(3,1);            % Velocity storage full
GE = u.grav;             % Gravity estimate

%init_x = [C; D; INT; VS; VSf; GE];
%% Simulating with Spine Lab Kalman Filter
x_1_1 = zeros(3,1);
sig_1_1 = params.sigAfferent * 10;

particle_states = randn( 6*nParticles, 1 ).*sigC(1,1);

full_states = zeros( 3, length(t) );
full_sigs = zeros( 3, length(t) );
full_particles = zeros( 6*nParticles, length(t) );

full_states(:,1) = x_1_1';
full_sigs(:,1) = diag( sig_1_1 );
full_particles(:,1) = particle_states';

x_prev = x_1_1;
sig_prev = sig_1_1;
particle_prev = particle_states;

figure(21); clf; hold on

for i=2:length( t )
    alpha = angAcc(:,i);
    [x_new, sig_new, particle_new] = SpineLabParticleStep( x_prev, sig_prev, particle_prev, alpha, params, dt, t(i) );
    
    full_states(:,i) = x_new';
    full_sigs(:,i) = diag( sig_new );
    full_particles(:,i) = particle_new';
    
    x_prev = x_new;
    sig_prev = sig_new;
    particle_prev = particle_new;
end

if mopo == 2
    % Positive and Negative Curve Fit
    full_rot = params.Tcan'*full_states;
    [m,i] = max( full_rot(2,:) );
    t1 = t(i:find( t > 19, 1, 'first' ) );
    y1 = full_rot(2,i:find( t > 19, 1, 'first' ) );
    
    [m,i] = min( full_rot(2,:) );
    t2 = t(i:end );
    y2 = full_rot(2,i:end );
    
    x0 = [y1(1), y2(1), 30];
    x1 = fminsearch(@(x) FitExponentialDecay( x, t1, t1, y1, y1 ), x0 )
    x2 = fminsearch(@(x) FitExponentialDecay( x, t2, t2, y2, y2 ), x0 )
end