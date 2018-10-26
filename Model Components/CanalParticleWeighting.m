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

function [fused_estimate, fused_sigma, afferent_weights] = CanalParticleWeighting( particle_afferents, all_means, all_sigs )
    
    nAfferents = size( particle_afferents, 2 );
    
    % Compute the new weights for the afferents
    afferent_weights = ones( nAfferents, 1 );
    for i=1:size( all_means, 1 )
        if sum(all_sigs(i,:) == 0) == 0
            afferent_weights = afferent_weights .* mvnpdf( particle_afferents', all_means(i,:), all_sigs(i,:) );
        end
        % Renormalize each time to prevent weights going to zero
        if sum( afferent_weights ) < 1e-10
            afferent_weights = ones( nAfferents, 1 ) / nAfferents;
        else
            afferent_weights = afferent_weights / sum( afferent_weights );
        end
    end
    
    % Weighted mean
    fused_estimate = ( particle_afferents * afferent_weights ); % / nAfferents;
    fused_sigma = zeros( size( fused_estimate ) );
    for i = 1:length(fused_estimate)
        fused_sigma(i) = sqrt( (particle_afferents(i,:) - fused_estimate(i)).^2 * afferent_weights );
    end
end