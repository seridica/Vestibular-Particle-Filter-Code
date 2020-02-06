%%%
% File: SpineLabDualAdaptationDeriv.m
% Author: Calvin Kuo
% Date: 10-12-2018
%
% This code uses two canal adaptation lines and fuses them via particle
% reweighting and a zero-mean gaussian prior on the afferent signal D.
% CAVEAT: At the moment, both long-term adaptation lines use the same canal
% input, but can easily add another canal line.
%
% States and their derivatives
%   - C = Canal Position (3x1)
%   - D1 = Canal Afferent First channel (3x1)
%   - D2 = Canal Afferent Second channel (3x1)
%   - INT = Endolymph position (3x1)
%   - VS = Velocity Leakage (3x1)
%   - VSf = Full Velocity Storage (3x1)
%   - GE = Gravitational Estimate (3x1)
%
% Other variables that are computed (you'll have to compute them again
% offline)
%   - Dobs = Bayesian combined D1 and D2 signals with zero-mean prior
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
%       : tn1 = Long term adaptation time constant for D1
%       : tn2 = Long term adaptation time constant for D2
%       : sigprior = Covariance on the afferent prior
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

function dx = SpineLabDualAdaptationDeriv( t, x, u, params )
    
    %% STATES
    C = x(1:3);
    D1 = x(4:6);
    D2 = x(7:9);
    INT = x(10:12);
    VS = x(13:15);
    VSf = x(16:18);
    GE = x(19:21);
    
    %% CANAL PROCESSING STEPS
    % Canal dynamics
    dC = CanalDynamicsDeriv( t, C, u, params );
    
    % Canal afferents long-term adaptation (set dD = dC if you want to
    % skip this)
    u.dC = dC;
    params.tn = params.tn1;
    dD1 = CanalAdaptationDeriv( t, D1, u, params );
    
    params.tn = params.tn2;
    dD2 = CanalAdaptationDeriv( t, D2, u, params );
    
    % Fuse adaptation lines
    [D, Sig, AW] = CanalParticleWeighting( [D1,D2], [0,0,0], params.sigprior );
    %D = mean( [D1, D2], 2 );
    
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
    dx = [dC; dD1; dD2; dINT; dVS; dVSf; dGE];
end