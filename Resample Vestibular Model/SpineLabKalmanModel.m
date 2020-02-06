%%%
% File: ParticleModel.m
% Author: Calvin Kuo
% Date: 09-10-2019
% Notes: Revised from the original. This is a Kalman-esque implementation
% of the resampling filter. By Kalman-esque, I simply mean that
% distributions are represented with Gaussians

%% Model parameters
% Simulation Parameters
params.tc = 4;
dt = 0.1;
Fs = 1/dt;

% Noise Sources
sigC = (1); % Noise on the input (Laurens and Angelaki)
sigA = (0.5); % Noise on the canal afferent
sigI = 0.05; % Internal noise
sigP = (5)^2; % Variance on the zero prior

params.sigAlpha = sigC;
params.sigAfferent = sigA;
params.sigState = sigI;
params.sigPrior = sigP;

%% Inputs
% Motion Profile
mopo = 2;
[t, angAcc] = MotionProfile( mopo, dt );
input = angAcc(2,:);

%% Simulating with Spine Lab Kalman Filter
figure(1); clf; hold on;
x_1_1 = zeros(6,1);
sig_1_1 = [sigC, sigC, sigA, sigA, 7^2, sigC]';

states = zeros( 6, length(t) );
sigs = zeros( 6, length(t) );

states(:,1) = x_1_1;
sigs(:,1) = sig_1_1;

x_prev = x_1_1;
sig_prev = sig_1_1;

for i=2:length( t )
    alpha = input(i);
    [x_new, sig_new] = SpineLabKalmanStep( x_prev, sig_prev, alpha, params, dt );
    
    states(:,i) = x_new;
    sigs(:,i) = sig_new;
    
    x_prev = x_new;
    sig_prev = sig_new;
end

figure(1); clf; hold on;
plot( t, states(1,:) );
plot( t, states(2,:) );
plot( t, states(3,:) );
plot( t, states(4,:) );
plot( t, states(5,:) );
plot( t, states(6,:)*5 );

figure(2); clf; hold on
plot( t, sigs(5,:) )

percep = states(5,:);
uncert = sigs(5,:);
t_onstart = find( t > 2, 1, 'first' );
t_onend = find( t > 60, 1, 'first' );

t_offstart = find( t > 62, 1, 'first' );
t_offend = find( t > 120, 1, 'first' );

t1 = t(t_onstart:t_onend);
t2 = t(t_offstart:t_offend);

y1 = states(6,t_onstart:t_onend);
y2 = states(6,t_offstart:t_offend);
x0 = [1, 1, 10];
x1 = fminsearch(@(x) FitExponentialDecay( x, t1, t1, y1, y1 ), x0 )
x2 = fminsearch(@(x) FitExponentialDecay( x, t2, t2, y2, y2 ), x0 )


% 
% if mopo == 2
%     % Positive and Negative Curve Fit
%     full_rot = params.Tcan'*states(13:15,:) + states(1:3,:);
%     [m,i] = max( full_rot(2,:) );
%     t1 = t(i:find( t > 45, 1, 'first' ) );
%     y1 = full_rot(2,i:find( t > 45, 1, 'first' ) );
%     
%     [m,i] = min( full_rot(2,:) );
%     t2 = t(i:end );
%     y2 = full_rot(2,i:end );
%     
%     x0 = [y1(1), y2(1), 30];
%     x1 = fminsearch(@(x) FitExponentialDecay( x, t1, t1, y1, y1 ), x0 )
%     x2 = fminsearch(@(x) FitExponentialDecay( x, t2, t2, y2, y2 ), x0 )
% elseif mopo == 3
%     % Positive and Negative Curve Fit
%     full_rot = params.Tcan'*states(13:15,:);
%     [m,i] = min( full_rot(2,:) );
%     t1 = t(i:end );
%     y1 = full_rot(2,i:end );
%     
%     x0 = [y1(1), y1(1), 30];
%     x1 = fminsearch(@(x) FitExponentialDecay( x, t1, t1, y1, y1 ), x0 )
% end