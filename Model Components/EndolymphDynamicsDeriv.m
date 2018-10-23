%%%
% File: EndolymphDynamicsDeriv.m
% Author: Calvin Kuo
% Date: 10-10-2018
% Notes: This code is a restructuring of the Laurens_NatureNeuroscience_M.m
% code originally written by PAF, CJD, NK, and JSB
%
% This code provides specifically the endolymph dynamics equation
% States (x):
%    - INT = (1:3) endolymph position
% Parameters (params):
%    - kv = Gain on the canal afferent
% Additional inputs (u):
%    - D = Canal afferent signal (either from canal dynamics of from slow
%    adaptation)
%    - dVSnoise = The current VSnoise derivative
% Model:
%    - dINT = kv*D + dVSnoise

function dx = EndolymphDynamicsDeriv( t, x, u, params )
    % Make sure we have three states (the three endolymph states)
    assert( size( x, 1 ) == 3 );
    
    % Extract parameters
    kv = params.kv;
    
    % Check to see what kind of input we have. If it is a time series, then
    % we should have a 4xn matrix input. If it is already a single angular
    % acceleration value, then we should have a 3x1 input
    
    % Case when we are passed the angular acceleration at this time
    % directly
    vec_D = u.D;
    if ( size(vec_D,1) == 3 )
        assert( size(vec_D,2) == 1 );
        D = vec_D;
    % Case when we are passed a time series of angular acceleration and
    % must extract the current angular acceleration through interpolation
    elseif ( size(vec_D,1) == 4 )
        tD = vec_D(1,:);
        vD = vec_D(2:4,:);
        D = interp1( tD, vD', t )';
    % Something else, through an assertion
    else
        assert(1);
    end
    
    % Case when we are passed the derivative of the endolymph noise at this
    % time directly
    vec_dEnoise = u.dVSnoise;
    if ( size(vec_dEnoise,1) == 3 )
        assert( size(vec_dEnoise,2) == 1 );
        dEnoise = vec_dEnoise;
    % Case when we are passed a time series of angular acceleration and
    % must extract the current angular acceleration through interpolation
    elseif ( size(vec_dEnoise,1) == 4 )
        tdEnoise = vec_dEnoise(1,:);
        vdEnoise = vec_dEnoise(2:4,:);
        dEnoise = interp1( tdEnoise, vdEnoise', t )';
    % Something else, through an assertion
    else
        assert(1);
    end
    
    % Model
    dx = kv*D + dEnoise;
end