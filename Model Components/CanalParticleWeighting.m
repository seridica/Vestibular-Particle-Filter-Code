%%%
% File: CanalParticleWeighting.m
% Author: Calvin Kuo
% Date: 10-12-2018
%
% This code fuses canal afferent "particles" assuming a zero-mean gaussian
% prior on the afferent state
% Inputs:
%     - Afferents (3xn) - n afferent particles
%     - Prior covariance (3x3) - Covariance on the prior
% Outputs:
%     - Fused afferent (3x1)

function [fused_estimate, afferent_weights] = CanalParticleWeighting( particle_afferents, prior_sig )
    
    nAfferents = size( particle_afferents, 2 );
    
    % Compute the new weights for the afferents
    afferent_weights = mvnpdf( particle_afferents', [0,0,0], prior_sig );
    if sum( afferent_weights ) < 1e-5
        afferent_weights = ones( nAfferents, 1 ) / nAfferents;
    else
        afferent_weights = afferent_weights / sum( afferent_weights );
    end
    
    % Weighted mean
    fused_estimate = ( particle_afferents * afferent_weights ); % / nAfferents;
    
end