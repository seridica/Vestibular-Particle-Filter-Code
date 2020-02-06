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

nState = 200;
nAff = 50;
params.tc = ones(nAff,1)*4;

% For dual line
%params.tn = [ones(nAff/2,1)*25; ones(nAff/2,1)*10000];
params.tn = [ones(nAff/2,1)*10000; ones(nAff/2,1)*10000];

dt = 0.1;
Fs = 1/dt;

% Nois Histories
sigC = 0.1; %1.7; %[3.6];%[0.1]*pi/180 / dt;             % Canal noise (Laurens and Angelaki)
sigQ = 1.7; %1.7^2 * dt; %6.4^2; %0.01; %[14.0];

params.sigAlpha = sigC;
params.sigAfferent = sigQ;
params.sigPrior = 10.0;
%params.sigPrior = 10000.0;
params.sigState = 0.1;

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
%% Initializing the variables
sx_dist = randn( nState, 1) * sigQ;
sy_dist = randn( nState, 1) * sigQ;
sz_dist = randn( nState, 1) * sigQ;

cx_dist = zeros( nAff,1);
cy_dist = zeros( nAff,1);
cz_dist = zeros( nAff,1);

ax_dist = randn( nAff, 1 ) * sigQ;
ay_dist = randn( nAff, 1 ) * sigQ;
az_dist = randn( nAff, 1 ) * sigQ;

state_dist = {sx_dist, sy_dist, sz_dist, cx_dist, cy_dist, cz_dist, ax_dist, ay_dist, az_dist};

figure(22); clf; hold on;
mean_percep = zeros( length(t), 3 );

for i=2:length( t )
    alpha = angAcc(:,i);
    new_state_dist = SpineLabResampleStep( state_dist, alpha, params, dt, t(i) );
    
    %plot( i, [min( new_state_dist{2} ), max( new_state_dist{2} )], 'rx' );
    %plot( i, mean( new_state_dist{2} ), 'bo', 'MarkerSize', 5, 'LineWidth', 2 );
    
    plot( t(i), [min( new_state_dist{2} ), max( new_state_dist{2} )], 'rx' );
    plot( t(i), mean( new_state_dist{2} ), 'bo', 'MarkerSize', 5, 'LineWidth', 2 );
    if ( mod(i,20) == 0 )
        plot( t(i), new_state_dist{8}, 'ys' );
    end
    plot( t(i), mean( new_state_dist{5} ), 'ko', 'MarkerSize', 5, 'LineWidth', 2 );
    
    state_dist = new_state_dist;
    mean_percep(i,:) = [mean( new_state_dist{1} ), mean( new_state_dist{2} ), mean( new_state_dist{3} )];
    
    if mod(i,100) == 0
        currTime = t(i)
    end
end

plot( [t(1), t(end)], [0,0], 'k-.', 'LineWidth', 5 );

if mopo == 2
    % Positive and Negative Curve Fit
    full_rot = params.Tcan'*mean_percep';
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
    full_rot = params.Tcan'*mean_percep';
    [m,i] = min( full_rot(2,:) );
    t1 = t(i:end );
    y1 = full_rot(2,i:end );
    
    x0 = [y1(1), y1(1), 30];
    x1 = fminsearch(@(x) FitExponentialDecay( x, t1, t1, y1, y1 ), x0 )
end