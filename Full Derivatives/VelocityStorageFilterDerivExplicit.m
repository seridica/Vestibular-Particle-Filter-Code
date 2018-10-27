%%%
% File: VelocityStorageFilterDerivExplicit.m
% Author: Calvin Kuo
% Date: 10-26-2018
%
% This code is similar to the Karmali et. al. filter
%
% States and their derivatives
%   - C = Canal Position (3x1)
%   - D1 = Canal Afferent First channel (3x1)
%   - D2 = Canal Afferent Second channel (3x1)
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

function [dx, new_store] = VelocityStorageFilterDerivExplicit( t, x, store, u, params )
    
    %% STATES
    nChannels = length( x ) / 3 / 8;
    dx = zeros( size(x) );
    
    %% CANAL PROCESSING STEPS
    part_sens = zeros( nChannels*2, 3 ); % Particles for the canals
    
    %% EXTERNAL INPUTS
    for i=1:nChannels
        
        %% CHANNELS FROM THE ACTUAL SENSORS
        ooff = (i-1)*12;
        C1 = x( ooff+1:ooff+3 );
        dC1 = CanalDynamicsDeriv( t, C1, u, params );
        dC1 = AddNoise( t, dC1, u.Qnoise );
        
        % Canal afferents long-term adaptation (set dD = dC if you want to
        % skip this)
        D1 = x( ooff+4:ooff+6 );
        u.dC = dC1;
        params.tn = params.tn1;
        dD1 = CanalAdaptationDeriv( t, D1, u, params );
        part_sens((i-1)*2+1,:) = AddNoise(t, D1, u.Cnoise );
        
        C2 = x( ooff+7:ooff+9 );
        dC2 = CanalDynamicsDeriv( t, C2, u, params );
        dC2 = AddNoise( t, dC2, u.Qnoise );
        
        % Canal afferents long-term adaptation (set dD = dC if you want to
        % skip this)
        D2 = x( ooff+10:ooff+12 );
        u.dC = dC2;
        params.tn = params.tn2;
        dD2 = CanalAdaptationDeriv( t, D2, u, params );
        part_sens(i*2,:) = AddNoise( t, D2, u.Cnoise );
        
        dx((i-1)*12+1:i*12) = [dC1; dD1; dC2; dD2];
        
    end
    
    %% INTERNAL MODEL
    prev_corr_particles = reshape( store(7:(6+nChannels*2*3)), nChannels*2, 3 );
    sens_sigma = diag( std( part_sens ) ).^2;
    prev_sigma = diag( std( prev_corr_particles ) ).^2 + u.Qnoise.^2;
    dt = 0.1;
    
    corrected_particles = zeros( size( prev_corr_particles ) );
    for i=1:nChannels
        ooff = nChannels*3*4+(i-1)*12;
        C1_internal = x(ooff+1:ooff+3);
        D1_internal = x(ooff+4:ooff+6);
        corrected_particles((i-1)*2+1,:) =  inv(prev_sigma) * ( sens_sigma*((1-dt/params.tn1)*D1_internal - dt/params.tc*C1_internal - prev_corr_particles((i-1)*2+1,:)' + u.Qnoise*randn(3,1)) + prev_sigma*(part_sens((i-1)*2+1,:)') );
        
        C2_internal = x(ooff+7:ooff+9);
        D2_internal = x(ooff+10:ooff+12);
        corrected_particles(i*2,:) = inv(prev_sigma) * ( sens_sigma*((1-dt/params.tn1)*D2_internal - dt/params.tc*C2_internal - prev_corr_particles(i*2,:)' + u.Qnoise*randn(3,1) ) + prev_sigma*(part_sens(i*2,:)') );
    end
    
    %[omega_est, sigma_est, ww] = CanalParticleWeighting( corrected_particles', [0,0,0], params.sigprior );
    omega_est = mean( corrected_particles )';
    sigma_est = std( corrected_particles )';
    
    for i=1:nChannels
        C1_internal = x(ooff+1:ooff+3);
        intmod.alpha = ( corrected_particles((i-1)*2+1,:) - prev_corr_particles((i-1)*2+1,:) )' / dt;
        intmod.Tcan = eye(3);
        dC1 = CanalDynamicsDeriv( t, C1_internal, intmod, params );
        
        intmod.dC = AddNoise( t, dC1, u.Qnoise );
        D1_internal = x(ooff+4:ooff+6);
        params.tn = params.tn1;
        dD1 = CanalAdaptationDeriv( t, D1_internal, intmod, params );
        
        C2_internal = x(ooff+7:ooff+9);
        intmod.alpha = ( corrected_particles(i*2,:) - prev_corr_particles(i*2,:) )' / dt;
        dC2 = CanalDynamicsDeriv( t, C2_internal, intmod, params );
        
        intmod.dC = AddNoise( t, dC2, u.Qnoise );
        D2_internal = x(ooff+10:ooff+12);
        params.tn = params.tn2;
        dD2 = CanalAdaptationDeriv( t, D2_internal, intmod, params );
        
        dx(nChannels*2*3+(i-1)*12+1:nChannels*2*3+i*12) = [dC1; dD1; dC2; dD2];
    end
    
    %% Derivative
    new_store = [omega_est; sigma_est; reshape( corrected_particles, nChannels*2*3, 1 )];
end