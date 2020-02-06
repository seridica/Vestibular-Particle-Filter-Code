%%%
% File: RaphanAndCohenDeriv.m
% Author: Calvin Kuo
% Date: 8-23-2019
% Notes: This code recreates the Raphan and Cohen model in ode45 form.
% There is a separate file that uses state-space block-diagram approach and
% should provide the same results.

function x_next = RaphanAndCohenDeriv( x_pre, alpha, params, dt )
    
    % Extract parameters
    tc = params.tc; % Canal time constant
    tvs = params.tvs; % Velocity storage time constant
    g = params.g; % Gain for leaky integrator
    
    % State Vector
    canal_pre = x_pre(1:3);
    vel_pre = x_pre(4:6);
    
    % Canal State Dynamics
    A = diag( (1-dt./tc) );
    B = eye(3) * dt;
    canal_next = A * canal_pre + B * alpha;
    
    % Velocity storage dynamics
    Av = diag( (1-dt./tvs) );
    Bv = eye(3) * dt;
    vel_next = Av * vel_pre + Bv * g * canal_next;
    
    x_next = [canal_next; vel_next];
end