%%%
% File: KarmaliParticleDeriv.m
% Author: Calvin Kuo
% Date: 02-06-2020
% Notes: Computes the particle derivative

function [x_next, est_next, K_next] = KarmaliParticleDeriv( x_pre, canal_in, K_in, params, dt, st )
    
    % Extract parameters
    N = params.N; % Number of particles
    R = params.R; % Afferent variance
    Q = params.Q; % Input variance
    tc = params.tc; % Primary time constant
    tc2 = params.tc2; % High pass filter time constant

    % Find the variance of the previous state
    Pt = cov( x_pre ) / dt / 10;
    
    % Afferent variance
    Rt = R(1);
    
    % State model matrices
    A = [0, 1; -1/(tc(1)*tc2(1)), -(1/tc(1) + 1/tc2(1))];
    B = [0; 1];
    C = [0, 1/tc2(1)];
    
    % If stabilizing
    if st == 0
        Kt = [0;5];
    else
        Kt = Pt * C' / Rt; % Current Kalman gain
        Kt(1) = 0; % Set the first value to zero
    end
    K_next = [Kt, K_in(:,1:(end-1))];
    K = mean( K_next, 2 ); % Average over past Kalman Gains
    
    % Propagate particles
    x_next = zeros( size( x_pre ) );
    est_next = zeros( N, 1 );
    for i=1:N
        internal_prev = x_pre(i,:)'; % Previous particle state
        internal_deriv = A * internal_prev + K*(canal_in(i) - C*internal_prev);
        internal_next = internal_deriv * dt + internal_prev;
        
        x_next(i,:) = internal_next';
        est_next(i) = K(2) * (canal_in(i) - C*internal_prev);
    end

end