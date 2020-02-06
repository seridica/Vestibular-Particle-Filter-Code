%%%
% File: AccelerationEstimate.m
% Author: Calvin Kuo
% Date: 10-12-2018
% Notes: This code is a restructuring of the Laurens_NatureNeuroscience_M.m
% code originally written by PAF, CJD, NK, and JSB
%
% This code provides specifically the equation for velocity estimation
% Additional inputs (u):
%    - GIA = Gravito-inertial acceleration
%    - GE = gravity estimate
% Model:
%    - AE = GIA - GE

function AE = AccelerationEstimate( t, u, params )
    
    % Check to see what kind of input we have. If it is a time series, then
    % we should have a 4xn matrix input. If it is already a single angular
    % acceleration value, then we should have a 3x1 input
    
    % Case when we are passed the gravito-inertial acceleration at this
    % time directly
    vec_GIA = u.GIA;
    if ( size(vec_GIA,1) == 3 )
        assert( size(vec_GIA,2) == 1 );
        curr_GIA = vec_GIA;
    % Case when we are passed a time series of gravito-inertial
    % acceleration and must extract the current gravito-inertial
    % acceleration through interpolation
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
    % Case when we are passed a time series of gravitational estimate and
    % must extract the current gravitational estimate through
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
    
    % Model
    AE = GIA - GE;
end