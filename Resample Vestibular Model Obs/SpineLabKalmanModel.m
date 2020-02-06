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
%       : go = Direct visual pathway gain (set to 0 to turn off)
%       : ko = Indirect visual pathway gain (set to 0 to turn off)
%   - Velocity Estimate
%       : gd = Gain on canal afferent

params.tc1 = 3;
params.tc2 = 7;

% For dual line
params.tn1 = 10000; %25; %25; %75.9; %75;
params.tn2 = 10000;

dt = 0.1;
Fs = 1/dt;

% Nois Histories
sigC = eye(3)*(1)^2; %1.7; %[3.6];%[0.1]*pi/180 / dt;             % Canal noise (Laurens and Angelaki)
sigQ = eye(3)*(1)^2; %1.7^2 * dt; %6.4^2; %0.01; %[14.0];

params.sigAlpha = sigC;
params.sigAfferent = sigQ;
params.sigPrior = eye(3)*10000.0^2;
params.sigState = eye(3)*10.0^2;

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

% Angular position (for now)
params.Tcan = [1 0 0;...
          0 1 0;...
          0 0 1];

% uaxi = [sqrt(2)/2, 0, -sqrt(2)/2];
% cang = acos( sqrt(3)/3 );
% params.Tcan = [ ...
%     cos(cang) + uaxi(1)^2*(1-cos(cang)), uaxi(1)*uaxi(2)*(1-cos(cang)) - uaxi(3)*sin(cang), uaxi(1)*uaxi(3)*(1-cos(cang)) + uaxi(2)*sin(cang); ...
%     uaxi(2)*uaxi(1)*(1-cos(cang)) + uaxi(3)*sin(cang), cos(cang) + uaxi(2)^2*(1-cos(cang)), uaxi(2)*uaxi(3)*(1-cos(cang)) - uaxi(1)*sin(cang); ...
%     uaxi(3)*uaxi(1)*(1-cos(cang)) - uaxi(2)*sin(cang), uaxi(3)*uaxi(2)*(1-cos(cang)) + uaxi(1)*sin(cang), cos(cang) + uaxi(3)^2*(1-cos(cang)) ];

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

%% Simulating with Spine Lab Kalman Filter
figure(22); clf; hold on;
x_1_1 = zeros(15,1);%[randn(12,1);zeros(3,1)];
sig_1_1 = eye(3) * 10^2; %params.sigAfferent*20;

states = zeros( 15, length(t) );
sigs = zeros( 3, length(t) );

states(:,1) = x_1_1;
sigs(:,1) = ones(3,1)*100;

x_prev = x_1_1;
sig_prev = sig_1_1;

plot( t(1), x_1_1(2), 'ko' );
plot( [t(1), t(1)], [x_1_1(2) - sig_1_1(2), x_1_1(2) + sig_1_1(2)], 'k' );

for i=2:length( t )
    alpha = angAcc(:,i);
    [x_new, sig_new] = SpineLabKalmanStep( x_prev, sig_prev, alpha, params, dt, t(i) );
    
    states(:,i) = x_new;
    sigs(:,i) = diag( sig_new );
    
    x_prev = x_new;
    sig_prev = sig_new;
end

if mopo == 2
    % Positive and Negative Curve Fit
    full_rot = params.Tcan'*states(13:15,:) + states(1:3,:);
    [m,i] = max( full_rot(2,:) );
    t1 = t(i:find( t > 45, 1, 'first' ) );
    y1 = full_rot(2,i:find( t > 45, 1, 'first' ) );
    
    [m,i] = min( full_rot(2,:) );
    t2 = t(i:end );
    y2 = full_rot(2,i:end );
    
    x0 = [y1(1), y2(1), 30];
    x1 = fminsearch(@(x) FitExponentialDecay( x, t1, t1, y1, y1 ), x0 )
    x2 = fminsearch(@(x) FitExponentialDecay( x, t2, t2, y2, y2 ), x0 )
elseif mopo == 3
    % Positive and Negative Curve Fit
    full_rot = params.Tcan'*states(13:15,:);
    [m,i] = min( full_rot(2,:) );
    t1 = t(i:end );
    y1 = full_rot(2,i:end );
    
    x0 = [y1(1), y1(1), 30];
    x1 = fminsearch(@(x) FitExponentialDecay( x, t1, t1, y1, y1 ), x0 )
end