%%%
% File: LaurensParticleDeriv.m
% Author: Calvin Kuo
% Date: 11-3-2019
% Notes: This code solves the derivative for Laurens and Droulez's particle
% filter.

function x_next = LaurensParticleDeriv( x_pre, V, params, dt )

  tc = params.tc; % Canal time constant
    eta = params.eta; % Canal noise standard deviation
    omv = params.omv; % Prior standard deviation on the head angular velocity
    N = params.N; % Number of particles
    
    % Size of x_pre indicates number of particles
    assert( size( x_pre, 1 ) == N );
    canal_pre = x_pre(:,1);
    omega_pre = x_pre(:,2);

    % Assertions for bookkeeping
    assert( length( eta ) == N );
    assert( length( omv ) == N );
    %assert( length( tc ) == N );
    
    % Generate noise vector
    canal_noise = normrnd( 0, eta );
    
    % Compute noisy current canal based on observation
    canal_curr = V + canal_noise;
    
    % Compute the current head velocity estimate based on inverse canal
    % dynamics model
    omega_curr = omega_pre - ( 1 - dt./tc ) .* canal_pre + canal_curr;
    
    x_next = LaurensReweight( [ canal_curr, omega_curr ], omv );
end

function x_resampled = LaurensReweight( x_raw, pdist )
    raw_weights = normpdf( x_raw(:,2), 0, pdist );
    norm_weights = raw_weights / sum( raw_weights );
    cum_weights = cumsum( norm_weights );
    
    x_resampled = zeros( size( x_raw ) );
    for i=1:length( x_raw )
        x_resampled(i,:) = x_raw( find( rand() <= cum_weights, 1, 'first' ), : );
    end
end