%%%
% File: LaurensVestibularModelDeriv.m
% Author: Calvin Kuo
% Date: 10-10-2018
% Notes: This code is a restructuring of the Laurens_NatureNeuroscience_M.m
% code originally written by PAF, CJD, NK, and JSB
%
% This code differs from the original in that it only provides the
% derivative of various states that are tracked in the model.
% This allows this code to be used in conjunction with standard matlab
% integrators.
%
% Note, the individual models for components are broken out into separate
% files so they can be used interchangeably and for other purposes (e.g.
% the canal model can be used to drive particles in a particle filter
% implementation, or can be doubled to make two separate lines).
%
% States and their derivatives
%   - C = Canal Position (3x1)
%   - D = Canal Afferent (3x1)
%   - INT = Endolymph position (3x1)
%   - VS = Velocity Leakage (3x1)
%   - VSf = Full Velocity Storage (3x1)
%   - GE = Gravitational Estimate (3x1)
%
% Other variables that are computed (you'll have to compute them again
% offline)
%   - rSL = Retinal Slip
%   - indVI = Indirect visual pathway
%   - dirVI = Direct visual pathway
%   - VE = Velocity estimate
%   - GIA = Gravito-inertial acceleration
%
% Inputs (all in ground frame)
%   - Angular velocity (constant 3x1 or time varying 4xn)
%   - Angular acceleration (constant 3x1 or time varying 4xn)
%   - Gravity (constant 3x1 or time varying 4xn)
%   - Linear acceleration (constant 3x1 or time varying 4xn)
%   - Canal Noise (constant 3x1 or time varying 4xn)
%   - Endolymph Noise Derivative (constant 3x1 or time varying 4xn)
%   - Visual Noise (constant 3x1 or time varying 4xn)
%   - Canal orientation (constant 3x3 or time varying 4x3xn)
%
% Parameters
%   - Canal Dynamics
%       : tc = Canal time constant
%       : tn = Long term adaptation time constant
%       : kv = Gain on canal afferent
%       : tvs = Velocity leakage time constant
%   - Gravity Estimation
%       : kf = Gain on gravity correction (set to 0 to turn off)
%       : ts = Gravity estimate time constant
%   - Visual Feedback
%       : go = Direct visual pathway gain (set to 0 to turn off)
%       : ko = Indirect visual pathway gain (set to 0 to turn off)
%   - Velocity Estimate
%       : gd = Gain on canal afferent

function dx = LaurensVestibularModelDeriv( t, x, u, params )
    
    %% STATES
    C = x(1:3);
    D = x(4:6);
    INT = x(7:9);
    VS = x(10:12);
    VSf = x(13:15);
    GE = x(16:18);
    
    %% CANAL PROCESSING STEPS
    % Canal dynamics
    dC = CanalDynamicsDeriv( t, C, u, params );
    
    % Canal afferents long-term adaptation (set dD = dC if you want to
    % skip this)
    u.dC = dC;
    dD = CanalAdaptationDeriv( t, D, u, params );
    
    % Endolymph motion
    u.D = D;
    dINT = EndolymphDynamicsDeriv( t, INT, u, params );
    
    %% VISION PROCESSING
    % Leaky Velocity Storage
    u.indVI = [0;0;0];
    u.GE = [1;0;0];
    u.GIA = [1;0;0];
    u.dINT = dINT;
    tempV = params.vswitch;
    tempG = params.gswitch;
    params.vswitch = 0;
    params.gswitch = 0;
    dVS = VelocityStorage( t, VS, u, params );
    params.vswitch = tempV;
    params.gswitch = tempG;
    
    % Retinal Slip and Visual Pathways
    u.VS = VS;
    [rSL, indVI, dirVI] = RetinalSlip( t, u, params );
    
    %% GRAVITY PROCESSING
    % Gravito-inertial acceleration
    GIA = GravitInertialAccel( t, u, params );
    
    % Velocity Estimate
    u.dirVI = dirVI;
    u.VSf = VSf;
    VE = VelocityEstimate( t, u, params );
    
    % Gravity estimate
    u.GIA = GIA;
    u.VE = VE;
    dGE = GravityEstimate( t, GE, u, params );
    
%     %% VELOCITY STORAGE
    % Full Velocity Storage
    u.indVI = indVI;
    u.GE = GE;
    dVSf = VelocityStorage( t, VSf, u, params );
%     
%     %% Derivative
    dx = [dC; dD; dINT; dVS; dVSf; dGE];
    %dx = [dC; dD; dINT; dVS; dVSf];
end