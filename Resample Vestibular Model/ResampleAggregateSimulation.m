%%%
% File: ResampleAggregateSimulation.m
% Author: Calvin Kuo
% Date: 11-18-2019
% Notes: This code runs simulations a few times for a specific settings and
% averages and fits the results.

function [t, output_estimate, output_zero, output_std] = ResampleAggregateSimulation( dur, nSims, debug, make_plots )
    
    % Simulation Parameters
    nState = 250;
    nAff = 50;
    dt = 0.1;
    Fs = 1/dt;
    
    % Model parameters
    params.tc = ones( nAff, 1 ) * 4;
    params.tn = ones( nAff, 1 ) * 10000;
    params.sigAlpha = 1.0; % Input variance
    params.sigAfferent = 1.0; % Canal afferent variance
    params.sigPrior = 25.0; % Zero prior for perception
    params.sigState = 0.01; % Process variance
    
    % Motion profile generation
    mopo = 2;
    [t, angAcc3] = MotionProfile( mopo, dt, dur );
    angAcc = angAcc3(2,:);
    
    % Setup data structures
    output_estimate = zeros(nSims, length(t));
    output_zero = zeros(nSims, length(t));
    output_std = zeros(nSims, length(t));
    
    % Simulate
    for i=1:nSims
        
        % Initializing the variables
        inter_dist = randn( nState, 1 ) * 2.0;
        canal_dist = zeros( nAff, 1 );
        aff_dist = randn( nAff, 1 ) * sqrt( params.sigAfferent );

        state_dist = {inter_dist, canal_dist, aff_dist};
        
        if make_plots == 1
            canal_state = zeros( 1, length(t) );
        end
        
        % Simulation run
        for j=1:length(t)
            alpha = angAcc(j);
            new_state_dist = SpineLabResampleStep( state_dist, alpha, params, dt, t(j) );
            output_estimate(i,j) = mean( new_state_dist{1} );
            output_zero(i,j) = mean( new_state_dist{4} );
            output_std(i,j) = std( new_state_dist{1} );
            
%             if mod(j,100) == 0
%                 simNum = i
%                 currTime = t(j)
%             end
           
            if make_plots == 1
                canal_state(j) = mean( new_state_dist{2} );
            end
            
            state_dist = new_state_dist;
        end
        
        if make_plots == 1
            figure(1); clf; hold on;
            plot( t, canal_state, 'k' );
            plot( t, output_estimate(i,:), 'b' );
            plot( t, output_estimate(i,:) + output_std(i,:), 'r' );
            plot( t, output_estimate(i,:) - output_std(i,:), 'r' );
            plot( t, output_zero(i,:), 'g' );
        end
        
        if debug == 1
            keyboard;
        end
    end
end