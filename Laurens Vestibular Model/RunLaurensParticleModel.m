%%%
% File: RunLaurensParticleModel.m
% Author: Calvin Kuo
% Date: 02-06-2020
% Notes: Function for running the Laurens Particle Model. Only returns the
% perceptual state. Assumes there is only one rotational acceleration axis
% as the input

function percep = RunLaurensParticleModel( t, angAcc, params, rep, seed )

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
    
    %% Simualte true canal dynamics
    Vfull = zeros(length(t), 1);
    for i=2:length(t)
        Vfull(i) = ( ( -1 ./ params.tc ) .* Vfull(i-1) - angAcc(i) ) .* dt + Vfull(i-1);
    end
    
    %% Laurens Particle Filter
    % Output perception storage
    percep = zeros(length(t), 1);
    
    % Initialize the state
    x_curr = zeros(params.N,2);
    
    for i=1:length(t)
        V = Vfull(i);
        
        % Next time step
        x_new = LaurensParticleDeriv( x_curr, V, params, dt );
        
        % Store
        percep(i) = mean( x_new(:,2) );
        x_curr = x_new;
        
        % Report status
        if rep == 1 && mod(i,Fs) < 1.0
            ['Laurens Model ', num2str(t(i))]
        end
    end
end