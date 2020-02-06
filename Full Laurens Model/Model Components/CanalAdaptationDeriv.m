%%%
% File: CanalAdaptationDeriv.m
% Author: Calvin Kuo
% Date: 10-10-2018
% Notes: This code is a restructuring of the Laurens_NatureNeuroscience_M.m
% code originally written by PAF, CJD, NK, and JSB
%
% This code provides specifically the slow adaptation dynamics of the canal
% afferents
% States (x):
%    - D = Canal afferent adaptation state
% Parameters (params):
%    - tn = the time constant for the long term adaptation
% Additional inputs (u):
%    - dC = Canal state (ideally from CanalDynamicsDeriv)
% Model:
%    - dD = -D/Tn + dC

function dx = CanalAdaptationDeriv( t, x, u, params )
    % Make sure we have three states (the three canal states)
    assert( size( x, 1 ) == 3 );
    
    % Extract parameters
    tn = params.tn;
    
    % Check to see what kind of input we have. If it is a time series, then
    % we should have a 4xn matrix input. If it is already a single angular
    % acceleration value, then we should have a 3x1 input
    
    % Case when we are passed the canal velococity at this time
    % directly
    vec_u = u.dC;
    if ( size(vec_u,1) == 3 )
        assert( size(vec_u,2) == 1 );
        dC = vec_u;
    % Case when we are passed a time series of canal velocities and
    % must extract the current canal velocity through interpolation
    elseif ( size(vec_u,1) == 4 )
        tdC = vec_u(1,:);
        vdC = vec_u(2:4,:);
        dC = interp1( tdC, vdC', t )';
    % Something else, through an assertion
    else
        assert(1);
    end
    
    % Model
    dx = -x/tn + dC;
end