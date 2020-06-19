%%%
% File: RunObserverModel.m
% Author: Calvin Kuo
% Date: 02-06-2020
% Notes: Function for running the Observer Model. Only returns the
% perceptual state. Assumes there is only one rotational acceleration axis
% as the input

function percep = RunObserverModel( t, angAcc, params, rep, seed )

    % Set the random number generator seed if requested
    if nargin == 5
        rng(seed);
    end
    
    if nargin < 4
        rep = 0;
    end
    
    % Get the time step and frequency based on the input t
    dt = ( t(end) - t(1) ) / length(t);
    Fs = 1/dt;
    
    % Make a 3-axis angular acceleration
    fullAngAcc = [angAcc'; zeros(2, length(t))];
    
    %% Raphan and Cohen Model
    % Output perception storage
    percep = zeros(length(t), 1);
    
    % Initialize the state
    x_prev = zeros(9,1);
    
    for i=2:length(t)
        % Angular acceleration input
        alpha = fullAngAcc(:,i);
        
        % Next time step
        x_new = ObserverDeriv( x_prev, alpha, params, dt );
        
        % Store
        percep(i) = x_new(7);
        x_prev = x_new;
        
        % Report status
        if rep == 1 && mod(i,Fs) < 1.0
            ['Observer Model ', num2str(t(i))]
        end
    end
end