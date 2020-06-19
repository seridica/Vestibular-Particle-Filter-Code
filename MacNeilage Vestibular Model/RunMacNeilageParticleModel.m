%%%
% File: RunMacNeilageParticleModel.m
% Author: Calvin Kuo
% Date: 02-06-2020
% Notes: Function for running the MacNeilage Particle Model. Only returns the
% perceptual state. Assumes there is only one rotational acceleration axis
% as the input

function percep = RunMacNeilageParticleModel( t, angAcc, params, rep, seed )

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
    
    %% MacNeilage Particle Filter
    % Output perception storage
    percep = zeros(length(t), 1);
    
    % Initialize the state
    x_curr = zeros(4,params.N);
    x_curr(4,:) = normrnd( 0, sqrt( params.R ) );
    
    for i=1:length(t)
        
        % Let state stabilize with a fixed Kalman Gain before releasing
        if ( i*dt<1.0 )
            st = 0;
        else
            st = 1;
        end
        
        % Next time step
        x_new = MacNeilageParticleDeriv( x_curr, angAcc(i), params, dt, st );
        
        % Store
        percep(i) = mean( x_new(2,:) );
        x_curr = x_new(1:4,:);
        
        % Report status
        if rep == 1 && mod(i,Fs) < 1.0
            ['Macneilage Model ', num2str(t(i))]
        end
    end
end