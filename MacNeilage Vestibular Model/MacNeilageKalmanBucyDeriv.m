%%%
% File: MacNeilageKalmanBucyDeriv.m
% Author: Calvin Kuo
% Date: 11-12-2019
% Notes: This code solves the derivative for MacNeilage's Kalman-Bucy 
% filter for just the semicircular canal and head orientations states
% This filter is loosely based on the Borah filter

function x_next = MacNeilageKalmanBucyDeriv( x_pre, V, params, dt )
    
    % Parameter extract
    Q = params.Q;
    R = params.R; % * dt; % / 0.001;
    tc = params.tc;
    beta = params.beta;
    
    % State and state variance extract
    canal_state_pre = x_pre(1);
    internal_state_pre = x_pre(2:4);
    est_var_pre = reshape( x_pre(5:13), 3, 3 );

    % State transition matrix
    A = [0, 1, 0; ...
         -beta^2, -2*beta, 0; ...
         0, 1, -tc];
     
    % Canal transition matrix
    Ac = -tc;
    Bc = 1;
     
    % Input matrix
    E = [0; 1; 0];
    
    % Observation matrix
    C = [0, 0, 1];
    
    % Power spectral density scaling for sensor noise
    %rho = 0.001;

    % Sensor derivative
    canal_state_deriv = Ac * canal_state_pre + Bc * V;
    
%     % Internal state variance derivative
%     internal_state_var_deriv = A * internal_state_var_pre + internal_state_var_pre * A' + E * Q * E';
%     
%     % Sensor variance
%     R = rho * diag( C * internal_state_var_pre * C' );
%     
    % Estimate Variance derivative
    est_var_deriv = A * est_var_pre + est_var_pre * A' + E * Q * E' - est_var_pre * C' * inv(R) * C * est_var_pre;
    %est_var_deriv = A * est_var_pre + est_var_pre * A' - est_var_pre * C' * inv(R) * C * est_var_pre;
    
    % Kalman Gain
    GK = est_var_pre * C' * inv(R);
        
%     if abs( V ) > 1e-1
%         keyboard;
%     end
    
    % Internal State derivative
    internal_state_deriv = A * internal_state_pre + GK * (canal_state_pre - C * internal_state_pre );
    
    % Update internal state
    internal_state_next = internal_state_pre + dt * internal_state_deriv;
%     
%     % Update internal state variance
%     internal_state_var_next = internal_state_var_pre + dt * internal_state_var_deriv;
    
    % Update estimate state variance
    est_var_next = est_var_pre + dt * est_var_deriv;
    
    % Update canal state
    canal_state_next = canal_state_pre + dt * canal_state_deriv;
    
    % Output
    x_next = [canal_state_next; internal_state_next; reshape( est_var_next, 9, 1 )];
    %keyboard;
end