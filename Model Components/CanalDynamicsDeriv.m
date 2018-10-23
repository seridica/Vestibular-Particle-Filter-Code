%%%
% File: CanalDynamicsDeriv.m
% Author: Calvin Kuo
% Date: 10-10-2018
% Notes: This code is a restructuring of the Laurens_NatureNeuroscience_M.m
% code originally written by PAF, CJD, NK, and JSB
%
% This code provides specifically the dynamics of the canals
% States (x):
%    - C = Canal afferent state
% Parameters (params):
%    - tc = the time constant for the canal
% Additional inputs (u):
%    - alpha = angular acceleration input
%    - noiseC = canal noise (this needs to be an input to avoid issues with
%    ode integration of random variables)
%    - Tcan = rotation matrix from lab frame to canal alignment
% Model:
%    - dC = -C/Tc + alpha + noise (optional)
%
% Note, to get rid of noise, simply set sigC param to sigC = zeros(3,3) or
% sigC = 0 (use matrix covariance if you want different canals to have
% different noise characteristics and cross-coupled errors).

function dx = CanalDynamicsDeriv( t, x, u, params )
    % Make sure we have three states (the three canal states)
    assert( size( x, 1 ) == 3 );
    
    % Extract parameters
    tc = params.tc;
    
    % Check to see what kind of input we have. If it is a time series, then
    % we should have a 4xn matrix input. If it is already a single angular
    % acceleration value, then we should have a 3x1 input
    
    % Case when we are passed the angular acceleration at this time
    % directly
    vec_alpha = u.alpha;
    if ( size(vec_alpha,1) == 3 )
        assert( size(vec_alpha,2) == 1 );
        curr_alpha = vec_alpha;
    % Case when we are passed a time series of angular acceleration and
    % must extract the current angular acceleration through interpolation
    elseif ( size(vec_alpha,1) == 4 )
        talpha = vec_alpha(1,:);
        alpha = vec_alpha(2:4,:);
        curr_alpha = interp1( talpha, alpha', t )';
    % Something else, through an assertion
    else
        assert(1);
    end
    
    % Case when we are passed the noise at this time
    % directly
    vec_noise = u.Cnoise;
    if ( size(vec_noise,1) == 3 )
        assert( size(vec_noise,2) == 1 );
        Cnoise = vec_noise;
    % Case when we are passed a time series of angular acceleration and
    % must extract the current angular acceleration through interpolation
    elseif ( size(vec_noise,1) == 4 )
        tnoise = vec_noise(1,:);
        all_noise = vec_noise(2:4,:);
        
        % Interp or floor
        Cnoise = interp1( tnoise, all_noise', t )';
    % Something else, through an assertion
    else
        assert(1);
    end
    
    % Case when we are passed the Tcan at this time
    % directly
    vec_Tcan = u.Tcan;
    if ( length( size(vec_Tcan) ) == 2 )
        Tcan = vec_Tcan;
    % Case when we are passed a time series of Tcan and
    % must extract the current Tcan through interpolation (tricky)
    % FOR NOW - CHOOSE THE CLOSEST TCAN NEED TO FIX
    elseif ( length( size(vec_Tcan) == 3 ) )
        Tcan_time = reshape( vec_Tcan(1,:,:), size( vec_Tcan, 3 ), 1 );
        tind = find( Tcan_time < t, 1, 'last' );
        Tcan = reshape( vec_Tcan(2:4,:,tind), 3, 3 );
    % Something else, through an assertion
    else
        assert(1);
    end
    
    % Model
    dx = -x/tc + Tcan*curr_alpha + Cnoise;
end