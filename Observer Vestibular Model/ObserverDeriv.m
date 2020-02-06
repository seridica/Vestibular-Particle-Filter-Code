%%%
% File: RaphanAndCohenDeriv.m
% Author: Calvin Kuo
% Date: 8-23-2019
% Notes: This code recreates the Raphan and Cohen model in ode45 form.
% There is a separate file that uses state-space block-diagram approach and
% should provide the same results.

function x_next = ObserverDeriv( x_pre, alpha, params, dt )
    
    % Extract parameters
    tc = params.tc; % Canal time constant
    K = params.K; % Gain for observer
    
    % State Vector
    canal_pre = x_pre(1:3);
    internal_canal_pre = x_pre(4:6);
    vel_pre = x_pre(7:9);
    
    % Canal State Dynamics
    A = diag( (1-dt./tc) );
    B = eye(3) * dt;
    canal_next = A * canal_pre + B * alpha;
    
    % Internal sensor model dynamics
    internal_canal_next = ( K * canal_next + A * internal_canal_pre - vel_pre ) / (K+1);
    
    % Observer
    vel_next = K * (canal_next - internal_canal_next);
    
    x_next = [canal_next; internal_canal_next; vel_next];
end