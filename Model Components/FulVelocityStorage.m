%%%
% File: FullVelocityStorage.m
% Author: Calvin Kuo
% Date: 10-10-2018
% Notes: This code is a restructuring of the Laurens_NatureNeuroscience_M.m
% code originally written by PAF, CJD, NK, and JSB
%
% This code provides specifically the dynamics of full velocity storage
% States (x):
%    - VSf = Full velocity storage
% Parameters (params):
%    - kf = Gain on the gravity estimate
% Additional inputs (u):
%    - GIA = gravito-inertial acceleration
%    - GE = Gravity estimate
%    - indVI = Indirect visual pathway
%    - dVS = leaky velocity storage
%    - Tcan = rotation matrix from lab frame to canal alignment
% Model:
%    - dVSf = dVS + indVI + kf*Tcan*( GIA/norm(GIA) x GE/norm(GE) )
%
% Note, to get rid of noise, simply set sigC param to sigC = zeros(3,3) or
% sigC = 0 (use matrix covariance if you want different canals to have
% different noise characteristics and cross-coupled errors).

function dx = FulVelocityStorage( t, x, u, params )
    % Make sure we have three states (the three canal states)
    assert( size( x, 1 ) == 3 );
    
    % Extract parameters
    kf = params.kf;
    
    % Check to see what kind of input we have. If it is a time series, then
    % we should have a 4xn matrix input. If it is already a single angular
    % acceleration value, then we should have a 3x1 input
    
    % Case when we are passed the gravito-inertial acceleration at this
    % time directly
    vec_GIA = u.GIA;
    if ( size(vec_GIA,1) == 3 )
        assert( size(vec_GIA,2) == 1 );
        curr_GIA = vec_GIA;
    % Case when we are passed a time series of angular acceleration and
    % must extract the current gravito-inertial acceleration through
    % interpolation
    elseif ( size(vec_GIA,1) == 4 )
        tGIA = vec_GIA(1,:);
        GIA = vec_GIA(2:4,:);
        curr_GIA = interp1( tGIA, GIA', t )';
    % Something else, through an assertion
    else
        assert(1);
    end
    
    % Case when we are passed the gravitational estimate at this time
    % directly
    vec_GE = u.GE;
    if ( size(vec_GE,1) == 3 )
        assert( size(vec_GE,2) == 1 );
        curr_GE = vec_GE;
    % Case when we are passed a time series of gravity estimate and
    % must extract the current gravitational acceleration through
    % interpolation
    elseif ( size(vec_GE,1) == 4 )
        tGE = vec_GE(1,:);
        all_GE = vec_GE(2:4,:);
        
        % Interp or floor
        curr_GE = interp1( tGE, all_GE', t )';
    % Something else, through an assertion
    else
        assert(1);
    end
    
    % Case when we are passed the leaky velocity storage derivative at this
    % time directly
    vec_dVS = u.dVS;
    if ( size(vec_dVS,1) == 3 )
        assert( size(vec_dVS,2) == 1 );
        curr_dVS = vec_dVS;
    % Case when we are passed a time series of leaky velocity storage
    % derivative and must extract the current leaky velocity storage
    % derivative through interpolation
    elseif ( size(vec_dVS,1) == 4 )
        tdVS = vec_dVS(1,:);
        all_dVS = vec_dVS(2:4,:);
        
        % Interp or floor
        curr_dVS = interp1( tdVS, all_dVS', t )';
    % Something else, through an assertion
    else
        assert(1);
    end
    
    % Case when we are passed the indirect visual pathway at this time
    % directly
    vec_indVI = u.indVI;
    if ( size(vec_indVI,1) == 3 )
        assert( size(vec_indVI,2) == 1 );
        curr_indVI = vec_indVI;
    % Case when we are passed a time series of indirect visual pathway and
    % must extract the current indirect visual pathway through interpolation
    elseif ( size(vec_indVI,1) == 4 )
        tindVI = vec_indVI(1,:);
        all_indVI = vec_indVI(2:4,:);
        
        % Interp or floor
        curr_indVI = interp1( tindVI, all_indVI', t )';
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
    dx = curr_dVS + curr_indVI + kf*Tcan*cross( curr_GIA / norm(curr_GIA), curr_GE / norm(curr_GE) );
end