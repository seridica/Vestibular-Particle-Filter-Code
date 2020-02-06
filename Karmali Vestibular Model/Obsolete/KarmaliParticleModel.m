%%%
% File: KarmaliParticleModel.m
% Author: Calvin Kuo
% Date: 8-26-2018
% Notes: Runs the Karmali Particle model

%% Model parameters
% Parameters
%   - Canal Dynamics
%       : tc = Canal time constant
%       : K = Observer gain

% From Karmali
dt = 1/1000;
%dt = 0.1;
Fs = 1/dt;

% Parameters from Karmali Paper
N = 33;
params.N = N;
params.tc = ones(1,N) * 5.7;
params.t2 = 0.005;
params.R = ones(1,N) * 2.8;
params.Q = ones(1,N) * 41;

%% Inputs
% Motion profile (angular acceleration)

% Motion Profile
mopo = 2;
[t, angAcc] = MotionProfile( mopo, dt );
angAcc = cumtrapz( angAcc(2,:) ) * dt;

%% Simulating with Observer Model - Deriv
x_prev = zeros(4,N);
%x_prev(1,:) = normrnd(0,sqrt(params.Q));
%x_prev(2,:) = normrnd(0,sqrt(params.Q));
%x_prev(3,:) = normrnd(0,sqrt(params.R));
%x_prev(4,:) = normrnd(0,sqrt(params.R));
full_states = zeros(8,length(t));
%full_states(1:4,1) = mean( x_prev' )';

for i=2:length( t )
    alpha = angAcc(i);
    if i > floor( 1.0/dt )
        st = 1;
    else
        st = 0;
    end
    x_new = KarmaliParticleDeriv( x_prev, alpha, params, dt, st );
    
    full_states(:,i) = mean( x_new' )';
    x_prev = x_new(1:4,:);
    if mod(i,10*Fs) == 0
        i*dt
    end
end

figure(1); clf; hold on;
plot( t, full_states(2,:), 'k-' );
plot( t, full_states(4,:), 'r-' );
plot( t, full_states(6,:), 'g-' );

figure(5); clf; hold on
plot( t, full_states(4,:) );

%% Simulating with Observer Model - TRANSFER FUNCTION
figure(2); clf;

canal_tf_x = tf([1], [1,1/params.tc(1)]);
canal_tf_y = tf([1], [1,1/params.tc(2)]);
canal_tf_z = tf([1], [1,1/params.tc(3)]);

canal_model_x = tf([1,0], [1,1/params.tc(1)]);
canal_model_y = tf([1,0], [1,1/params.tc(2)]);
canal_model_z = tf([1,0], [1,1/params.tc(3)]);

K = 3;

[y,t,x] = lsim([canal_tf_y; canal_model_y*(K*canal_tf_y) / (1+K*canal_model_y); (K*canal_tf_y) / (1+K*canal_model_y)], angAcc(2,:), t);
clf; hold on;
plot( t, y(:,1), 'k-' );
plot( t, y(:,2), 'r-' );
plot( t, y(:,3), 'g-' )
xlim([-50,150])