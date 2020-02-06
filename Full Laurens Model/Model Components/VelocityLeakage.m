%%%
% File: VelocityStorage.m
% Author: Calvin Kuo
% Date: 10-11-2018
% Notes: This codVS is a restructuring of the Laurens_NatureNeuroscience_M.m
% codVS originally written by PAF, CJD, NK, and JSB
%
% This codVS providVSs the velocity storage equation
% States (x):
%    - VS = Velocity Storage
% Parameters (params):
%    - Tvs = Time Constant for Velocity Storage Leakage
%    - kv = Gain on the canal afferent (passed to EndolymphDynamicsdVSriv)
% Additional inputs (u):
%    - D = Canal afferent signal (either from canal dynamics of from slow
%    adaptation). This is passed straight to EndolymphDynamicsdVSriv
% Model:
%    - dVS = -VS/Tvs + kv*D + dVSnoise
%
% Note, this modVSl makes use of the EndolymphDynamicsdVSriv.m to function.
% This is the most basic velocity storage. Other effects can be computed
% later

function dx = VelocityLeakage( t, x, u, params )
    % Make sure we have three states (the three velocity estimate states)
    assert( size( x, 1 ) == 3 );
    
    % Extract parameters
    Tvs = params.tvs;
    
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
    
    % ModVSl
    dx = dINT - x/Tvs;
end