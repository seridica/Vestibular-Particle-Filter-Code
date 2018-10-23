%%%
% File: GravitoInertialAccel.m
% Author: Calvin Kuo
% Date: 10-12-2018
% Notes: This code is a restructuring of the Laurens_NatureNeuroscience_M.m
% code originally written by PAF, CJD, NK, and JSB
%
% This code specifically provides calculation for GIA
% Additional inputs (u):
%    - grav = gravity vector (lab frame)
%    - acc = linear acceleration vector (lab frame)
% Model:
%    - GIA = grav - acc

function GIA = GravitInertialAccel( t, u, params )
    
    % Check to see what kind of input we have. If it is a time series, then
    % we should have a 4xn matrix input. If it is already a single angular
    % acceleration value, then we should have a 3x1 input
    
    % Case when we are passed the gravity acceleration at this time
    % directly
    vec_grav = u.grav;
    if ( size(vec_grav,1) == 3 )
        assert( size(vec_grav,2) == 1 );
        curr_grav = vec_grav;
    % Case when we are passed a time series of gravity acceleration and
    % must extract the current angular acceleration through interpolation
    elseif ( size(vec_grav,1) == 4 )
        tgrav = vec_grav(1,:);
        grav = vec_grav(2:4,:);
        curr_grav = interp1( tgrav, grav', t )';
    % Something else, through an assertion
    else
        assert(1);
    end
    
    % Case when we are passed the linear acceleration at this time
    % directly
    vec_acc = u.acc;
    if ( size(vec_acc,1) == 3 )
        assert( size(vec_acc,2) == 1 );
        curr_acc = vec_acc;
    % Case when we are passed a time series of angular acceleration and
    % must extract the current angular acceleration through interpolation
    elseif ( size(vec_acc,1) == 4 )
        tacc = vec_acc(1,:);
        all_acc = vec_acc(2:4,:);
        
        % Interp or floor
        curr_acc = interp1( tacc, all_acc', t )';
    % Something else, through an assertion
    else
        assert(1);
    end
    
    % Model
    GIA = curr_grav - curr_acc;
end