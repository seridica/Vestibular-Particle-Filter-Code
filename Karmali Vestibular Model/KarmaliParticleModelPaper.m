%%%
% File: KarmaliParticleModelPaper.m
% Author: Calvin Kuo
% Date: 8-27-2018
% Notes: Runs the Karmali Particle model exactly implemented from their
% paper

%% Model parameters
% Parameters
%   - Canal Dynamics
%       : tc = Canal time constant
%       : K = Observer gain

% From Karmali
dt = 1/8000;
%dt = 0.01;
Fs = 1/dt;

% Parameters from Karmali Paper
%N = 158;
N = 33;
%tc = normrnd( 5.7, 1.0, N, 1 ); %ones(N,1) * 5.7;
tc = ones(N,1) * 5.7;
tc2 = ones(N,1) * 0.005;
%R = [ones(N,1) * 13];% ones(158,1) * 2.8];
R = [ones(N,1) * 2.8];
%Q = [ones(N,1) * 190.0]; %[ones(33,1) * 41];% ones(158,1) * 190];
Q = [ones(N,1)*41];

%% Inputs
% Motion profile (angular acceleration)

% Motion Profile
mopo = 5;
[t, angAcc] = MotionProfile( mopo, dt );
angVel = cumtrapz( angAcc(2,:) )*dt;
tlen = length(t);

%% Simulate Canals with Noise
canal_outs = zeros( N, tlen );
canal_temps = zeros( N, tlen );
[b,a] = butter(1, 400/(Fs/2));

for i=1:N
    i
    
    % Generate noise
    input_noise_raw = normrnd( 0, sqrt(Q(i)), [1,tlen] );
    input_noise_filt = filter( b, a, input_noise_raw );
    input_noise_filt = input_noise_filt * (sqrt(Q(i)) / std( input_noise_filt ));
    %input_noise_filt = input_noise_filt * (Q(i) / std( input_noise_filt ));
    clear input_noise_raw;
    
    proc_noise_raw = normrnd( 0, sqrt(R(i)), [1,tlen] );
    proc_noise_filt = filter( b, a, proc_noise_raw );
    proc_noise_filt = proc_noise_filt * (sqrt(R(i)) / std( proc_noise_filt ));
    %proc_noise_filt = proc_noise_filt * (R(i) / std( proc_noise_filt ));
    clear proc_noise_raw;
    
    %canal_curr = normrnd( 0, sqrt(Q(i)), [2,1] );
    canal_curr = [0;0];
    canal_outs(i,1) = proc_noise_filt(1);
    for j=2:tlen
        A = [0, 1; -1/(tc(i)*tc2(i)), -(1/tc(i) + 1/tc2(i))];
        B = [0; 1];
        C = [0, 1/tc2(i)];
        canal_deriv = A * canal_curr + B * ( angVel(j) + input_noise_filt(j) );
        canal_curr = canal_deriv * dt + canal_curr;
        %canal_curr = ( A*dt + eye(2) ) * canal_curr + B * dt * ( angVel(j) + input_noise_filt(j) );
        canal_temps(i,j) = C*canal_curr;
        canal_outs(i,j) = C*canal_curr + proc_noise_filt(j);
    end
    clear proc_noise_filt
    clear input_noise_filt
    
end

%canal_outs = canal_outs * (4.0 / mean( std( canal_outs(i,1:Fs) ) ) );

figure(1); clf; hold on;
plot( t, mean( canal_outs ) );

mean( std( canal_outs(:,1:8000) ) ) % This should be 4.0 for regular afferents
std( std( canal_outs(:,1:800) ) ) % This should be 1.0 for regular afferents

%% Simulating with Particle Observer Model
%internal_curr = normrnd( 0, 1e-3, [N, 2] );
internal_curr = zeros( N, 2 );
internal_states_1 = zeros( N, tlen );
internal_states_2 = zeros( N, tlen );
internal_est = zeros( N, tlen );
internal_outs = zeros( N, tlen );
Ks = zeros( 4, tlen );
for i=2:tlen
    Pt = cov( internal_curr ) / dt / (Fs * 1.1 / 1000 + 1);
    %Pt = [0, 0; 0, var( internal_curr(:,2) ) / dt];
    %Pt = std( internal_states_2 )^2 / dt; % Note, only using the second state
    %Rt = std( canal_outs(:,i) )^2;
    Rt = R(1); %R(1);
    
    if (mod(i,1/dt) == 0)
        i*dt
    end

    A = [0, 1; -1/(tc(1)*tc2(1)), -(1/tc(1) + 1/tc2(1))];
    B = [0; 1];
    C = [0, 1/tc2(1)];

    % Compute Kalman Gain
    if i < floor( 1.0 / dt )
        K = [0,5]';
        Kt = Pt * C' / Rt;
        %Ks(1:2,i) = Kt;
        Ks(1:2,i) = K;
        Ks(3:4,i) = K;
    else
        Kt = Pt * C' / Rt;
        %Kt = Pt * C' / R(j);
        %Kt = Pt * C' / Q(j);
        Ks(1:2,i) = Kt;
        K = mean( Ks(1:2,(i-floor(0.05/dt)):i)' )';
        K(1) = 0;
        %K = Kt;
        Ks(3:4,i) = K;
    end
    %K = [0,3]';
    
    %invM = inv( eye(2) + K*C*dt );
    %invM = inv( eye(2) + K*C );
        
    % Go through particles
    for j=1:N
        
        internal_prev = internal_curr(j,:)';
        %internal_curr = [internal_states(j); internal_states_2(j)];
        %vel_pre = K*( canal_outs(j,i-1) - C*internal_curr );
        %internal_curr = ( A*dt + eye(2) ) * internal_curr + B * dt * vel_pre;
        
        %internal_this = invM * ( (dt*A + eye(2)) * internal_this + dt*K*canal_outs(j,i) );
        
        %internal_this = (dt*A + eye(2)) * internal_this + dt*K*( canal_outs(j,i-1) - C*internal_this);
        %internal_next = (dt*A + eye(2) - K*C*dt) * internal_prev + dt*K*canal_outs(j,i-1);
        
        internal_deriv = A * internal_prev + K*(canal_outs(j,i-1) - C*internal_prev);
        internal_next = internal_deriv*dt + internal_prev;
        internal_outs(j,i) = K(2)*(canal_outs(j,i-1) - C*internal_prev);%C*internal_curr;
        
        internal_curr(j,1) = internal_next(1);
        internal_curr(j,2) = internal_next(2);
        
        internal_est(j,i) = C * internal_next; 
    end
    internal_states_1(:,i) = internal_curr(:,1);
    internal_states_2(:,i) = internal_curr(:,2);
end

plot( t, mean( internal_outs ), 'g' );
%mean( std( internal_outs ) ) % This should be 0.87 for regular afferents
std( mean( internal_outs(:,1:8000) ) ) % This is the correct threshold criteria

figure(5); clf; hold on
plot( t, Ks(2,:) );
plot( t, Ks(4,:) );

%% Simulating with Observer Model - TRANSFER FUNCTION
figure(2); clf;

canal_tf_x = tf([1], [1,1/tc(1)]);
canal_tf_y = tf([1], [1,1/tc(2)]);
canal_tf_z = tf([1], [1,1/tc(3)]);

canal_model_x = tf([1,0], [1,1/tc(1)]);
canal_model_y = tf([1,0], [1,1/tc(2)]);
canal_model_z = tf([1,0], [1,1/tc(3)]);

K = 3;
[y,t,x] = lsim([canal_tf_y; canal_model_y*(K*canal_tf_y) / (1+K*canal_model_y); (K*canal_tf_y) / (1+K*canal_model_y)], angAcc(2,:), t);
clf; hold on;
plot( t, y(:,1), 'k-' );
plot( t, y(:,2), 'r-' );
plot( t, y(:,3), 'g-' )
xlim([-50,150])