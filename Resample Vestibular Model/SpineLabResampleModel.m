%%%
% File: SpineLabResampleModel.m
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

if ~exist('plt','var')
    plt = 0;
end

nState = 250;
nAff = 50;
params.tc = ones(nAff,1)*4;
%params.tc = randn([nAff,1])*1.0 + 5;

% For dual line
params.tn = ones(nAff,1)*10000;

dt = 0.1;
Fs = 1/dt;

% Nois Histories
sigC = 1.0; %1.7; %[3.6];%[0.1]*pi/180 / dt;             % Canal noise (Laurens and Angelaki)
sigQ = 1.0; %1.7^2 * dt; %6.4^2; %0.01; %[14.0];

params.sigAlpha = sigC;
params.sigAfferent = sigQ;
params.sigPrior = 25.0;
params.sigState = 0.01;

% Profile duration
if ~exist('dur','var')
    dur = 40;
end

%% Inputs

% Motion Profile
mopo = 2;
[t, angAcc] = MotionProfile( mopo, dt, dur );

%% Initializing the variables
inter_dist = randn( nState, 1) * 2.0;
canal_dist = zeros( nAff,1);
aff_dist = randn( nAff, 1 ) * sqrt( sigQ );

state_dist = {inter_dist, canal_dist, aff_dist};

%figure(22); clf; hold on;
mean_percep = zeros( length(t), 1 );
zero_percep = zeros( length(t), 1 );
std_percep = zeros( length(t), 1 );

plot(t, cumtrapz( angAcc(2,:) )*0.1, 'Color', [0, 0.4470, 0.7410], 'LineWidth', 5);
for i=2:length( t )
    alpha = angAcc(2,i);
    new_state_dist = SpineLabResampleStep( state_dist, alpha, params, dt, t(i) );
    
    if plt && ( i/10 > 20-5 && i/10 < 20+dur+65 )
        %plot( i, [min( new_state_dist{2} ), max( new_state_dist{2} )], 'rx' );
        %plot( i, mean( new_state_dist{2} ), 'bo', 'MarkerSize', 5, 'LineWidth', 2 );
    
        %plot( t(i), mean( new_state_dist{1} ) + [-std( new_state_dist{1} ), std( new_state_dist{1} )], 'Color', [0.5,0.5,0.5], 'Marker', 'o' );
        plot( t(i), mean( new_state_dist{1} ), 'Color', [0.1294, 0.6667, 0.1294], 'Marker', 'o', 'MarkerSize', 5, 'LineWidth', 2 );
        if ( mod(i,5) == 0 )
           cur_dist = new_state_dist{1};
           plot( t(i), cur_dist(1:5:250), 'Marker', 'x', 'LineStyle', 'none', 'MarkerSize', 5, 'LineWidth', 0.5, 'Color', [0.8, 0.8, 0.8, 0.2] );
        end
        plot( t(i), mean( new_state_dist{2} ), 'Color', [0.8500, 0.3250, 0.0980], 'Marker', 'o', 'MarkerSize', 5, 'LineWidth', 2 );
    end
    
    state_dist = new_state_dist;
    mean_percep(i,:) = mean( new_state_dist{1} );
    zero_percep(i,:) = mean( new_state_dist{4} );
    std_percep(i,:) = std( new_state_dist{1} );
   
    if mod(i,100) == 0
        currTime = t(i)
    end
end

%plot( [t(1), t(end)], [0,0], 'k-.', 'LineWidth', 5 );

if mopo == 2
    % Positive and Negative Curve Fit
    i_on = find( t > 2.0, 1, 'first' );
    t_on = find( t > dur, 1, 'first' );
    t1 = t(i_on:t_on);
    ymean_on = mean_percep(i_on:t_on)';
    ystd_on = std_percep(i_on:t_on)' - mean( std_percep((end-10*Fs):end));
    yzero_on = zero_percep(i_on:t_on)';
    
    i_off = find( t > dur + 2, 1, 'first' );
    t_off = find( t > dur + 60, 1, 'first' );
    t2 = t(i_off:t_off );
    ymean_off = mean_percep(i_off:t_off)';
    ystd_off = std_percep(i_off:t_off)' - mean( std_percep((end-10*Fs):end));
    yzero_off = zero_percep(i_off:t_off)';
    
    xmean_0 = [ymean_on(1), ymean_off(1), 10];
    xmean_on = fminsearch(@(x) FitExponentialDecay( x, t1, t1, ymean_on, ymean_on ), xmean_0 )
    xmean_off = fminsearch(@(x) FitExponentialDecay( x, t2, t2, ymean_off, ymean_off ), xmean_0 )
    
    
    xstd_0 = [ystd_on(1), ystd_off(1), 10];
    xstd_on = fminsearch(@(x) FitExponentialDecay( x, t1, t1, ystd_on, ystd_on ), xstd_0 )
    xstd_off = fminsearch(@(x) FitExponentialDecay( x, t2, t2, ystd_off, ystd_off ), xstd_0 )
    
    
    xzero_0 = [yzero_on(1), yzero_off(1), 10];
    xzero_on = fminsearch(@(x) FitExponentialDecay( x, t1, t1, yzero_on, yzero_on ), xzero_0 )
    xzero_off = fminsearch(@(x) FitExponentialDecay( x, t2, t2, yzero_off, yzero_off ), xzero_0 )
   
    
    figure(5); clf; hold on
    plot( t, mean_percep, 'r' );
    plot( t1, xmean_on(1) * exp( -(t1-t1(1)) / xmean_on(3) ), 'r-.' );
    plot( t2, xmean_off(1) * exp( -(t2-t2(1)) / xmean_off(3) ), 'r-.' );
    
    plot( t, zero_percep, 'g' );
    plot( t1, xzero_on(1) * exp( -(t1-t1(1)) / xzero_on(3) ), 'g-.' );
    plot( t2, xzero_off(1) * exp( -(t2-t2(1)) / xzero_off(3) ), 'g-.' );
    
    plot( t, std_percep, 'b' );
    plot( t1, xstd_on(1) * exp( -(t1-t1(1)) / xstd_on(3) ) + mean( std_percep((end-10*Fs):end)), 'b-.' );
    plot( t2, xstd_off(1) * exp( -(t2-t2(1)) / xstd_off(3) ) + mean( std_percep((end-10*Fs):end)), 'b-.' );
end