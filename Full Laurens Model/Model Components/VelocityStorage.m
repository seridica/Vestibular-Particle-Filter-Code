%%%
% File: VelocityStorage.m
% Author: Calvin Kuo
% Date: 10-10-2018
% Notes: This code is a restructuring of the Laurens_NatureNeuroscience_M.m
% code originally written by PAF, CJD, NK, and JSB
%
% This code provides specifically the dynamics of full velocity storage
% States (x):
%    - VSf = Velocity storage
% Parameters (params):
%    - kf = Gain on the gravity estimate
%    - gswitch = include gravity (1 = yes, 0 = no)
%    - vswitch = incclude vision (1 = yes, 0 = no)
%    - Tvs = Time Constant for Velocity Storage Leakage
% Additional inputs (u):
%    - GIA = gravito-inertial acceleration
%    - GE = Gravity estimate
%    - indVI = Indirect visual pathway
%    - Tcan = rotation matrix from lab frame to canal alignment
%    - D = Canal afferent signal (either from canal dynamics of from slow
%    adaptation). This is passed straight to EndolymphDynamicsdVSriv
%    - dINT = change in endolymph position
% Model:
%    - dVS = -VS/Tvs + dINT + indVI + kf*Tcan*( GIA/norm(GIA) x GE/norm(GE) )

function dx = VelocityStorage( t, x, u, params )
    % Make sure we have three states (the three canal states)
    assert( size( x, 1 ) == 3 );
    
    % Extract parameters
    Tvs = params.tvs;
    kf = params.kf;
    gswitch = params.gswitch;
    vswitch = params.vswitch;
    
    % Check to see what kind of input we have. If it is a time series, then
    % we should have a 4xn matrix input. If it is already a single angular
    % acceleration value, then we should have a 3x1 input
    
    % Case when we are passed the endolymph velocity at this time
    % directly
    vec_u = u.dINT;
    if ( size(vec_u,1) == 3 )
        assert( size(vec_u,2) == 1 );
        dINT = vec_u;
    % Case when we are passed a time series of endolymph velocity and
    % must extract the current endolymph position through interpolation
    elseif ( size(vec_u,1) == 4 )
        tdINT = vec_u(1,:);
        vdINT = vec_u(2:4,:);
        dINT = interp1( tdINT, vdINT', t )';
    % Something else, through an assertion
    else
        assert(1);
    end
    
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
    dx = -x/Tvs + dINT + vswitch*curr_indVI + gswitch*kf*Tcan*cross( curr_GIA / norm(curr_GIA), curr_GE / norm(curr_GE) );
end