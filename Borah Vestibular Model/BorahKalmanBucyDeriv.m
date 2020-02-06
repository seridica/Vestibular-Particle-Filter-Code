%%%
% File: BorahKalmanBucyDeriv.m
% Author: Calvin Kuo
% Date: 11-12-2019
% Notes: This code solves the derivative for Borah's Kalman-Bucy filter for
% just the semicircular canal and head orientations tates

function x_next = BorahKalmanBucyDeriv( x_pre, V, params, dt )
    
    % Parameter extract
    Q = params.Q;
    
    % State and state variance extract
    internal_state_pre = x_pre(1:4);
    canal_state_pre = x_pre(5:8);
    internal_state_var_pre = reshape( x_pre(9:24), 4, 4 );
    est_var_pre = reshape( x_pre(25:40), 4, 4 );

    % State transition matrix
    A = [0, 1, 0, 0; ...
         -14400, -240, 0, 0; ...
         0, 0, 0, 1; ...
         0, 1, -0.00333, -0.1333];
     
    % Input matrix
    E = [0; 1; 0; 0];
    
    % Observation matrix
    C = [0, 0.5724, -0.0019, 57.17];
    
    % Power spectral density scaling for sensor noise
    rho = 0.001;

    % Sensor derivative
    canal_state_deriv = A * canal_state_pre + E * V;
%     
%     if abs( V ) > 1e-1
%         keyboard;
%     end
    
    % Sensor observation
    canal_obs = C * canal_state_pre;
    
    % Internal state variance derivative
    internal_state_var_deriv = A * internal_state_var_pre + internal_state_var_pre * A' + E * Q * E';
    
    % Sensor variance
    R = rho * diag( C * internal_state_var_pre * C' );
    
    % Estimate Variance derivative
    est_var_deriv = A * est_var_pre + est_var_pre * A' + E * Q * E' - est_var_pre * C' * inv(R) * C * est_var_pre;
    
    % Kalman Gain
    GK = est_var_pre * C' * inv(R);
    
    % Internal State derivative
    internal_state_deriv = A * internal_state_pre + GK * (canal_obs - C * internal_state_pre );
    
    % Update internal state
    internal_state_next = internal_state_pre + dt * internal_state_deriv;
    
    % Update internal state variance
    internal_state_var_next = internal_state_var_pre + dt * internal_state_var_deriv;
    
    % Update estimate state variance
    est_var_next = est_var_pre + dt * est_var_deriv;
    
    % Update canal state
    canal_state_next = canal_state_pre + dt * canal_state_deriv;
    
    % Output
    x_next = [internal_state_next; canal_state_next; reshape( internal_state_var_next, 16, 1 ); reshape( est_var_next, 16, 1 )];
    %keyboard;
end