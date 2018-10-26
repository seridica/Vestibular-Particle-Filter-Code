%%%
% File: VelocityStorageFilterDerivExplicit.m
% Author: Calvin Kuo
% Date: 10-23-2018
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

function [dx, new_store] = KarmaliDerivExplicit( t, x, store, u, params )
    
    %% STATES
    nChannels = length( x ) / 3 / 4;
    dx = zeros( size(x) );
    
    %% STORAGE
    prev_omega = store(1:3);
    prev_sigma = store(4:6);
    
    %% CANAL PROCESSING STEPS
    part_sens = zeros( nChannels*2, 3 );
    
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
    prev_corr_particles = reshape( store(10:(9+nChannels*2*3)), nChannels*2, 3 );
    prev_particles = reshape( store((10+nChannels*2*3):end), nChannels*2, 3 );

    Kgain = eye(3)*3; %diag( std( prev_particles ) ) * inv( u.Qnoise );
    corrected_particles = zeros( size( prev_corr_particles ) );
    for i=1:nChannels*2
        corrected_particles(i,:) = ( randn(1,3) * u.Qnoise + prev_corr_particles(i,:) + part_sens(i,:) - prev_particles(i,:)*(1-0.1/params.tc) ) * Kgain * inv( Kgain + eye(3) );
    end

    % Propagate internal model
    new_input = corrected_particles - prev_corr_particles;
    part_in = new_input + prev_particles*(1-0.1/params.tc);

    %[omega_est, sigma_est, ww] = CanalParticleWeighting( corrected_particles', [0,0,0], params.sigprior );
    omega_est = mean( corrected_particles )';
    sigma_est = std( corrected_particles )';
    
    %% Derivative
    new_store = [omega_est; sigma_est; diag(Kgain); reshape( corrected_particles, nChannels*2*3, 1 ); reshape( part_in, nChannels*2*3, 1 )];
end