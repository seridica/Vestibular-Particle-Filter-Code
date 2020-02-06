%%%
% File: MacNeilageParticleDeriv.m
% Author: Calvin Kuo
% Date: 11-15-2019
% Notes: This code solves the derivative for a particle filter that is
% based off of the MacNeilage Kalman-Bucy filter. Only the semicircular
% canal and head orientation states are included

function x_next = MacNeilageParticleDeriv( x_pre, V, params, dt, st )

    % Parameter extract
    Q = params.Q;
    R = params.R;
    tc = params.tc;
    beta = params.beta;
    N = params.N;
    
    % Extract relevant states (canal model and internal state estimate)
    canal_state_pre = x_pre(1,:);
    internal_state_pre = x_pre(2:4,:);
    
    % State Transition Matrix for Canal
    Ac = -tc;
    Bc = 1;
    
    % Generate noise
    proc_noise = normrnd( 0, sqrt(Q) );
    meas_noise = normrnd( 0, sqrt(R) );
    
    % Canal Propagate with noise
    noisy_input = V; % + proc_noise;
    canal_state_deriv = Ac .* canal_state_pre + Bc .* noisy_input;
    canal_state_next = canal_state_deriv * dt + canal_state_pre;
    noisy_canal = canal_state_next + meas_noise;
    
    % Estimate the covariance on the internal states
    Pt = cov( internal_state_pre' ) / dt;
    
    % Generate transitions for each particle (supports having slightly
    % different state transitions for each particle)
    internal_state_next = zeros( size( internal_state_pre ) );
    estimate_state = zeros( 1, N );
    for i=1:N
        A = [0, 1, 0; ...
             -beta(i)^2, -2*beta(i), 0; ...
             0, 1, -tc(i)];
        E = [0; 1; 0];
        C = [0, 0, 1];
        
        if st == 0
            K = [30;450;30];
        else
            K = Pt * C' / R(i);
        end
        
        % Force unobserved terms to 0
        % K(1:2) = 0;
        
        internal_state_deriv = A * internal_state_pre(:,i) + K * ( noisy_canal(i) - C * internal_state_pre(:,i) ) + E * proc_noise(i);
        internal_state_next(:,i) = internal_state_deriv * dt + internal_state_pre(:,i);
        estimate_state(i) = C * internal_state_next(:,i);
    end
    
    % Output
    x_next = [canal_state_next; internal_state_next; estimate_state];
end