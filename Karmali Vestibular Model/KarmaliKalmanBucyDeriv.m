%%%
% File: KarmaliKalmanBucyDeriv.m
% Author: Calvin Kuo
% Date: 11-13-2019
% Notes: This code implements Karmali's Kalman-Bucy filter

function x_next = KarmaliKalmanBucyDeriv( x_pre, V, params, dt )
    
    % Extract parameters
    R = params.R; % afferent variance
    Q = params.Q; % input variance
    tc = params.tc; % State transition parameter
    t2 = params.t2; % State transition parameter
    
    % Setup model matrices
    A = [0, 1; -1/(tc*t2), -( 1/tc + 1/t2 )];
    B = [0; 1];
    C = [0, 1/t2];
    
    % Extract previous states and variances
    prev_canal = x_pre(1:2);
    prev_state = x_pre(3:4);
    prev_var = reshape( x_pre(5:8), 2, 2 );
    
    % Canal derivatives (add noise?)
    canal_deriv = A*prev_canal + B*V;
    
    % Propagate canal
    next_canal = prev_canal + canal_deriv * dt;

    
    % Implicit Canal Integration
    %next_canal = inv( eye(2) - A * dt ) * (prev_canal + B*V*dt);

    % Compute Derivatives
    iR = inv(R);
    var_deriv = A*prev_var + prev_var*A' + B*Q*B' - prev_var*C'*iR*C*prev_var;
    next_var = prev_var + var_deriv * dt;
    
    % Implicit State Integration
    K = prev_var*C'*iR;
    state_deriv = A*prev_state + K*(C * prev_canal - C * prev_state);
    %next_state = inv(eye(2)-A*dt+K*C*dt) * (K*C*next_canal*dt + prev_state);
    
    % Propagate Derivatives
    next_state = prev_state + state_deriv * dt;
    next_obs = K * ( C * prev_canal - C * prev_state );
    
    % Fill posterior
    x_next = [next_canal; next_state; reshape( next_var, 4, 1 ); next_obs; K];
end