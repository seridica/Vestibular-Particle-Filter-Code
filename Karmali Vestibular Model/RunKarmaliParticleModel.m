%%%
% File: RunKarmaliParticleModel.m
% Author: Calvin Kuo
% Date: 02-06-2020
% Notes: Function for running the Karmali Particle Model. Only returns the
% perceptual state. Assumes there is only one rotational acceleration axis
% as the input

function percep = RunKarmaliParticleModel( t, angAcc, params, rep, seed )

    % Set the random number generator seed if requested
    if nargin == 5
        rng(seed);
    end
    
    if nargin < 4
        rep = 0;
    end
    
    % Get the time step and frequency based on the input t
    dt = ( t(end) - t(1) ) / length(t);
    Fs = 1/dt;
    
    % Get the angular velocity as the Karmali model uses the velocity as
    % input
    angVel = cumtrapz( angAcc )*dt;
    
    %% Simulate noisy canal dynamics
    % Pull necessary parameters
    N = params.N;
    Q = params.Q;
    R = params.R;
    tc = params.tc;
    tc2 = params.tc2;
    tlen = length(t);
    
    canal_outs = zeros( N, tlen );
    [b,a] = butter(1, 400/(Fs/2));

    for i=1:N
        % Generate input noise (filtered at 400Hz)
        input_noise_raw = normrnd( 0, sqrt(Q(i)), [1,tlen] );
        input_noise_filt = filter( b, a, input_noise_raw );
        input_noise_filt = input_noise_filt * (sqrt(Q(i)) / std( input_noise_filt ));
        clear input_noise_raw;

        % Generate process noise (filtered at 400Hz)
        proc_noise_raw = normrnd( 0, sqrt(R(i)), [1,tlen] );
        proc_noise_filt = filter( b, a, proc_noise_raw );
        proc_noise_filt = proc_noise_filt * (sqrt(R(i)) / std( proc_noise_filt ));
        clear proc_noise_raw;

        canal_curr = [0;0];
        canal_outs(i,1) = proc_noise_filt(1);
        for j=2:tlen
            A = [0, 1; -1/(tc(i)*tc2(i)), -(1/tc(i) + 1/tc2(i))];
            B = [0; 1];
            C = [0, 1/tc2(i)];
            canal_deriv = A * canal_curr + B * ( angVel(j) + input_noise_filt(j) );
            canal_curr = canal_deriv * dt + canal_curr;
            canal_outs(i,j) = C*canal_curr + proc_noise_filt(j);
        end
        clear proc_noise_filt
        clear input_noise_filt
    end
    
    %% Karmali Particle Filter
    % Storage for the Kalman gain
    K_store = [zeros( 1, floor( 0.05/dt ) ); ones(1, floor(0.05/dt)) * 5];
    
    % Output perception storage
    percep = zeros(length(t), 1);
    
    % Initialize the state
    x_curr = zeros(params.N,2);
    
    for i=2:length(t)
        % Get the current canal input
        canal_in = canal_outs(:,i-1);
        
        % Stabilize
        if i < floor( 1.0 / dt )
            st = 0;
        else
            st = 1;
        end
        
        % Next time step
        [x_new, est_out, K_new] = KarmaliParticleDeriv( x_curr, canal_in, K_store, params, dt, st );
        
        % Store
        percep(i) = mean( est_out );
        x_curr = x_new;
        K_store = K_new;
        
        % Report status
        if rep == 1 && mod(i,Fs) < 1.0
            ['Karmali Model ', num2str(t(i))]
        end
    end
end