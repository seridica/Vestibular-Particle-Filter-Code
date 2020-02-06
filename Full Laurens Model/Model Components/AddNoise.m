%%%
% File: AddNoise.m
% Author: Calvin Kuo
% Date: 10-23-2018
% Notes: This code is a restructuring of the Laurens_NatureNeuroscience_M.m
% code originally written by PAF, CJD, NK, and JSB
%
% This code specifically adds noise to some state at time t
% Model:
%   x_noise(t) = x_clean(t) + noise(t)
%
% Note, to get rid of noise, simply set sigC param to sigC = zeros(3,3) or
% sigC = zeros(3,1) (use matrix covariance if you want different canals to
% have different noise characteristics and cross-coupled errors).

function x_noise = AddNoise( t, x_clean, noise )
    % Find out how many states we have
    nStates = size( x_clean, 1 );
    
    % Determine how noise was given
    % Case where noise was given directly
    if ( size(noise,1) == nStates )
        
        % Noise was given as a 3-vector
        if ( size(noise,2) == 1 )
            add_noise = noise;
        % Noise was given as a covariance matrix - compute a random noise
        % now. WARNING: If using a varaible time step integrator, using
        % this option is dangerous and could lead to some slow integration
        % behavior
        else
            assert( size(noise,2) == nStates );
            add_noise = noise * randn(nStates,1);
        end
    % Case when we are passed a time series of noise values and
    % must extract the current noise through interpolation
    else
        assert( size(noise,1) == nStates+1 )
        tnoise = noise(1,:);
        all_noise = noise(2:4,:);
        
        % Interp or floor
        add_noise = interp1( tnoise, all_noise', t )';
    end
    
    % Add noise
    x_noise = x_clean + add_noise;
end