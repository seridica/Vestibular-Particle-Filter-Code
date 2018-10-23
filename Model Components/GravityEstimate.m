%%%
% File: GravityEstimate.m
% Author: Calvin Kuo
% Date: 10-11-2018
% Notes: This code is a restructuring of the Laurens_NatureNeuroscience_M.m
% code originally written by PAF, CJD, NK, and JSB
%
% This code provides the gravity estimate equation
% States (x):
%    - GE = Gravity estimate
% Parameters (params):
%    - Ts = Time Constant for translational acceleration decay
% Additional inputs (u):
%    - VE = Velocity estimate
%    - GIA = Gravito-inertial acceleration
%    - Tcan = Rotation to canal frame
% Model:
%    - dGE = GE x Tcan'*VE + (GIA - GE)/Ts

function dx = GravityEstimate( t, x, u, params )
    % Make sure we have three states (the three gravity estimate states)
    assert( size( x, 1 ) == 3 );
    
    % Extract parameters
    Ts = params.ts;
    
    % Case when we are passed the velocity estimate at this time
    % directly
    vec_VE = u.VE;
    if ( size(vec_VE,1) == 3 )
        assert( size(vec_VE,2) == 1 );
        VE = vec_VE;
    % Case when we are passed a time series of velocity estimates and
    % must extract the current velocity estimate through interpolation
    elseif ( size(vec_VE,1) == 4 )
        tVE = vec_VE(1,:);
        vVE = vec_VE(2:4,:);
        VE = interp1( tVE, vVE', t )';
    % Something else, through an assertion
    else
        assert(1);
    end
    
    % Case when we are passed the gravito-inertial acceleration at this
    % time directly
    vec_GIA = u.GIA;
    if ( size(vec_GIA,1) == 3 )
        assert( size(vec_GIA,2) == 1 );
        GIA = vec_GIA;
    % Case when we are passed a time series of gravito-inertial 
    % acceleration and must extract the current velocity estimate through
    % interpolation
    elseif ( size(vec_GIA,1) == 4 )
        tGIA = vec_GIA(1,:);
        vGIA = vec_GIA(2:4,:);
        GIA = interp1( tGIA, vGIA', t )';
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
    dx = (cross(x,Tcan'*VE)) - x/Ts + GIA/Ts;
end