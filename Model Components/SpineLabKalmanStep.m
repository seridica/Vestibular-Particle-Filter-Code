%%%
% File: SpineLabKalmanStep.m
% Author: Calvin Kuo
% Date: 11-02-2018

function [x_t_t, sig_t_t] = SpineLabKalmanStep( x_t1_t1, sig_t1_t1, alpha, params, dt )

    % Extract parameters
    tc = params.tc;
    tn1 = params.tn1;
    tn2 = params.tn2;
    
    sigAlpha = params.sigAlpha;
    sigAfferent = params.sigAfferent;
    sigPrior = params.sigPrior;

    % State Vector - Canal State (3), Long-Term State (3), Internal Model
    % State (3)
    canal_long1_t1 = x_t1_t1(1:6);
    canal_long2_t1 = x_t1_t1(7:12);
    int_t1_t1 = x_t1_t1(13:15);
    
    % Canal update and long-term state update
    A1 = zeros(6,6);
    A1(1:3,1:3) = eye(3)*(1-dt/tc);
    A1(4:6,1:3) = eye(3)*(-dt/tc);
    A1(4:6,4:6) = eye(3)*(1-dt/tn1);
    
    B1 = [eye(3); eye(3)] * dt;
    
    C1 = [eye(3)*dt, zeros(3,3); eye(3)*dt, eye(3)];
    
    canal_long1_t = A1 * canal_long1_t1 + B1 * params.Tcan * alpha; % + C1 * [sqrt(sigAlpha) * randn(3,1); sqrt(sigAfferent) * randn(3,1)];
    
    A2 = A1;
    A2(4:6,4:6) = eye(3)*(1-dt/tn2);
    
    canal_long2_t = A2 * canal_long2_t1 + B1 * alpha; % + C1 * [sqrt(sigAlpha) * randn(3,1); sqrt(sigAfferent) * randn(3,1)];
    
    % Internal state model, can update this as needed
    Aint = eye(3);
    int_t_t1 = Aint*int_t1_t1;
    sig_t_t1 = Aint'*sig_t1_t1*Aint + diag( eye(3)*abs( int_t_t1 - ( canal_long1_t(4:6) + canal_long2_t(4:6) ) / 2 ).^2 ) / 1000;
    %sig_t_t1 = diag( eye(3)*abs( int_t_t1 - ( canal_long1_t(4:6) + canal_long2_t(4:6) ) / 2 ).^2 )
    
    % Fusion
    sig_sensor = sigAlpha*dt + sigAfferent;
    denom = inv( sig_t_t1*sig_sensor*sig_sensor + sig_sensor*sig_sensor*sigPrior + 2*sig_sensor*sigPrior*sig_t_t1 );
    
    sig_t_t = sig_sensor * sig_sensor * sigPrior * sig_t_t1 * denom;
    int_t_t = denom * ( sig_sensor*sigPrior*sig_t_t1*canal_long1_t(4:6) + ...
                sig_sensor*sigPrior*sig_t_t1*canal_long2_t(4:6) + ...
                sig_sensor*sig_sensor*sigPrior*int_t_t1 );
            
    x_t_t = [ canal_long1_t; canal_long2_t; int_t_t ];
end