%%%
% File: SpineLabKalmanStep.m
% Author: Calvin Kuo
% Date: 9-10-2019

function [x_t_t, sig_t_t] = SpineLabKalmanStep( x_t1_t1, sig_t1_t1, alpha, params, dt )
    %[x_t_t, sig_t_t] = SpineLabKalmanStep1( x_t1_t1, sig_t1_t1, alpha, params, dt );
    [x_t_t, sig_t_t] = SpineLabKalmanStep2( x_t1_t1, sig_t1_t1, alpha, params, dt );
end

%% Original developed on 9-10
function [x_t_t, sig_t_t] = SpineLabKalmanStep1( x_t1_t1, sig_t1_t1, alpha, params, dt )
    % Get variance parameters
    sigC = params.sigAlpha;
    sigA = params.sigAfferent;
    sigI = params.sigState;
    sigP = params.sigPrior;
    
    % Get time constants
    tc = params.tc;
    
    % Compute canal dynamics
    canal_prev = x_t1_t1(2);
    canal_prev_sig = sig_t1_t1(2);
    A = 1 - dt/tc;
    B = dt;
    canal_curr = A * canal_prev + B * alpha;
    canal_curr_sig = A * canal_prev_sig * A' + B * sigC * B';
    
    % Input noise on the afferent
    canal_noise = canal_curr; % Afferent noise assumed to be zero mean
    canal_noise_sig = canal_curr_sig + sigA;
    
    % Resample
    state_prev_sig = sig_t1_t1(5);
    state_prior_sig = state_prev_sig + sigI;
    
    state_prev = x_t1_t1(5);
    state_prior = state_prev;
	
    noise_vec = [ 5000*canal_noise_sig; state_prior_sig ];
    noise_vec = noise_vec / sum( noise_vec );
    state_curr = noise_vec(1) * state_prior + noise_vec(2) * canal_noise;
    state_curr_sig = ( noise_vec(1) * ( state_prior_sig + state_prior^2 ) + ...
                     noise_vec(2) * ( canal_noise_sig + canal_noise^2 ) - ...
                     state_curr^2 );
    
    x_t_t = [alpha, canal_curr, canal_noise, state_prior, state_curr];
    sig_t_t = [sigC, canal_curr_sig, canal_noise_sig, state_prior_sig, state_curr_sig];
end

%% Developed 11-17
% Zero prior is applied on the state estimate and does not affect particle
% propagation

function [x_t_t, sig_t_t] = SpineLabKalmanStep2( x_t1_t1, sig_t1_t1, alpha, params, dt )
    % Get variance parameters
    sigC = params.sigAlpha;
    sigA = params.sigAfferent;
    sigI = params.sigState;
    sigP = params.sigPrior;
    
    % Get time constants
    tc = params.tc;
    
    % Compute canal dynamics
    canal_prev = x_t1_t1(2);
    canal_prev_sig = sig_t1_t1(2);
    A = 1 - dt/tc;
    B = dt;
    canal_curr = A * canal_prev + B * alpha;
    canal_curr_sig = A * canal_prev_sig * A' + B * sigC * B';
    
    % Input noise on the afferent
    canal_noise = canal_curr; % Afferent noise assumed to be zero mean
    canal_noise_sig = canal_curr_sig + sigA;
    
    % Propagate internal estimate
    state_prev_sig = sig_t1_t1(5);
    state_prior_sig = state_prev_sig + sigI;
    
    state_prev = x_t1_t1(5);
    state_prior = state_prev;
	
    % Resample
    noise_vec = [ 5000*canal_noise_sig; state_prior_sig ];
    noise_vec = noise_vec / sum( noise_vec );
    state_curr = noise_vec(1) * state_prior + noise_vec(2) * canal_noise;
    state_curr_sig = ( noise_vec(1) * ( state_prior_sig + state_prior^2 ) + ...
                     noise_vec(2) * ( canal_noise_sig + canal_noise^2 ) - ...
                     state_curr^2 );
    
    % Apply the zero prior with Bayesian
    zeroSigSum = state_curr_sig + sigP;
    state_est = sigP / zeroSigSum * state_curr + state_curr_sig / zeroSigSum * 0;
    state_est_sig = sigP * state_curr_sig / zeroSigSum;
    
    x_t_t = [alpha, canal_curr, canal_noise, state_prior, state_curr, state_est];
    sig_t_t = [sigC, canal_curr_sig, canal_noise_sig, state_prior_sig, state_curr_sig, state_est_sig];
end