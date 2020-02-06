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
    
    fused_estimate = zeros( size( all_means, 2 ), 1 );
    fused_sigma = zeros( size( all_means, 2 ), 1 );
    
    % Compute the new weights for the afferents
    afferent_weights = ones( nAfferents, size( all_means, 2 ) );
        
    % Treat components separately
    for i=1:size( all_means, 2 )

        for j=1:size( all_means, 1 )
            if sum(all_sigs(j,i) == 0) == 0
                afferent_weights(:,i) = afferent_weights(:,i) .* normpdf( particle_afferents(i,:)', all_means(j,i), sqrt( all_sigs(j,i) ) );
            end
            % Renormalize each time to prevent weights going to zero
            if sum( afferent_weights(:,i) ) < 1e-16
                if j == 1
                    afferent_weights(:,i) = ones( nAfferents, 1 ) / nAfferents;
                else
                    sum( afferent_weights(:,i) );
                    afferent_weights(:,i) = zeros( nAfferents, 1 );
                    [m, pind] = min( ( particle_afferents(i,j) - all_means(j,i) ).^2 );
                    afferent_weights(pind,i) = 1;
                    %afferent_weights(:,i) = ones( nAfferents, 1 ) / nAfferents;
                end
            else
                afferent_weights(:,i) = afferent_weights(:,i) / sum( afferent_weights(:,i) );
            end
        end

        % Weighted mean
        fused_estimate(i) = ( particle_afferents(i,:) * afferent_weights(:,i) ); % / nAfferents;
        %fused_sigma(i) = 0.2; %max( [( (particle_afferents(i,:) - fused_estimate(i) ).^2 * afferent_weights ), 0.2] );
        fused_sigma(i) = (particle_afferents(i,:) - fused_estimate(i) ).^2 * afferent_weights(:,i);
    end
end