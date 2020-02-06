%%%
% File: SimpleLaurensVestibularModelDerivExplicit.m
% Author: Calvin Kuo
% Date: 10-23-2018
% Notes: This code is a restructuring of the Laurens_NatureNeuroscience_M.m
% code originally written by PAF, CJD, NK, and JSB
%
% This code differs from the original in that it only provides the
% derivative of various states that are tracked in the model.
% This allows this code to be used in conjunction with standard matlab
% integrators. This particular derivative function only looks at the
% vestibular -> velocity estimate path (no vision, no gravity)
%
% States and their derivatives
%   - C = Canal Position (3x1)
%   - D = Canal Afferent (3x1)
%   - INT = Endolymph position (3x1)
%   - VS = Velocity Leakage (3x1)
%
% Inputs (all in ground frame)
%   - Angular acceleration (constant 3x1 or time varying 4xn)
%   - Canal Noise (constant 3x1 or time varying 4xn)
%   - Endolymph Noise Derivative (constant 3x1 or time varying 4xn)
%   - Canal orientation (constant 3x3 or time varying 4x3xn)
%
% Parameters
%   - Canal Dynamics
%       : tc = Canal time constant
%       : tn = Long term adaptation time constant
%       : kv = Gain on canal afferent
%       : tvs = Velocity leakage time constant
%   - Velocity Estimate
%       : gd = Gain on canal afferent

function [dx, new_store] = SimpleLaurensVestibularModelDerivExplicit( t, x, store, u, params )
    
    %% STATES
    C = x(1:3);
    D = x(4:6);
    INT = x(7:9);
    VS = x(10:12);
    
    %% CANAL PROCESSING STEPS
    % Canal dynamics
    dC = CanalDynamicsDeriv( t, C, u, params );
    dC = AddNoise( t, dC, u.Cnoise );
    
    % Canal afferents long-term adaptation (set dD = dC if you want to
    % skip this)
    u.dC = dC; % Add noise to the canal post integration
    dD = CanalAdaptationDeriv( t, D, u, params );
    
    % Endolymph motion
    u.D = D;
    dINT = EndolymphDynamicsDeriv( t, INT, u, params );
    dINT = AddNoise( t, dINT, u.VSnoise ); % Add noise on integration
    
    %% VISION PROCESSING
    % Leaky Velocity Storage
    u.indVI = [0;0;0];
    u.GE = [1;0;0];
    u.GIA = [1;0;0];
    u.dINT = dINT;
    params.vswitch = 0;
    params.gswitch = 0;
    dVS = VelocityStorage( t, VS, u, params );
    
    %% Derivative
    dx = [dC; dD; dINT; dVS];
    new_store = [];
end