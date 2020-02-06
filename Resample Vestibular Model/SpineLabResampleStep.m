%%%
% File: SpineLabResampleStep.m
% Author: Calvin Kuo
% Date: 02-22-2019
%
% Notes - This performs a step of the resample model for vestibular
% integration. Very similar to SpineLabParticleStep

function new_distribution = SpineLabResampleStep( state_distributions, alpha, params, dt, tnow )
    
    % Extract parameters
    tc = params.tc;
    tn = params.tn;
    
    sigAlpha = params.sigAlpha;
    sigAfferent = params.sigAfferent;
    sigPrior = params.sigPrior;
    sigState = params.sigState;
    
    % In the state distribution, the first set is the internal state
    % distribution and the second is the canal distribution
    internal_distribution = state_distributions{1};
    canal_distribution = state_distributions{2};
    afferent_distribution = state_distributions{3};
    
    % Update canals
    act_alpha = repmat( alpha, length( canal_distribution ), 1 ) + randn( [length( canal_distribution ), 1] ) * sqrt( sigAlpha );
    new_canal_distribution = canal_distribution .* ( 1 - dt./tc ) + act_alpha * dt;
    new_afferent_distribution = afferent_distribution .* ( 1 - dt./tn ) - canal_distribution .* dt./tc + act_alpha(:,1) * dt;
    
    % Resample afferents with the zero prior
    prior_afferent_distribution = new_afferent_distribution + randn( length( afferent_distribution ), 1) * sqrt(sigAfferent);
    %noisy_afferent_distribution = new_afferent_distribution + randn( length( afferent_distribution ), 1) * sigAfferent;
    %prior_afferent_distribution = ZeroPriorResample( noisy_afferent_distribution, sigPrior );
    %prior_afferent_distribution = new_afferent_distribution;
    
    % Propagate canals
%     if (tnow > 60 )
%         keyboard;
%     end
    resample_distributions = {internal_distribution, prior_afferent_distribution};
    new_internal_distribution = ResampleDistribution( resample_distributions, length(internal_distribution), 4, [500, 1] ) + sqrt(sigState) * randn( length( internal_distribution), 1 );
    %new_internal_distribution = ResampleDistribution( resample_distributions, length(internal_distribution), 4, [20, 1] );
    %new_internal_distribution = ZeroPriorResample( new_internal_distribution, sigPrior );
    zero_internal_distribution = ZeroPriorResample( new_internal_distribution, sqrt(sigPrior) );
    
    new_distribution = {new_internal_distribution, ...
                        new_canal_distribution, ...
                        new_afferent_distribution, ...
                        zero_internal_distribution};
    %plot( tnow, mean( prior_afferent_distribution ), 'co', 'MarkerSize', 5, 'LineWidth', 2  );
    %plot( tnow, mean( prior_afferent_distribution ) + std( prior_afferent_distribution ), 'cd', 'MarkerSize', 5, 'LineWidth', 2  );
    %plot( tnow, mean( prior_afferent_distribution ) - std( prior_afferent_distribution ), 'cd', 'MarkerSize', 5, 'LineWidth', 2  );
end

function outDist = ZeroPriorResample( inDistribution, zeroSigma )
    probVec = normpdf( inDistribution, 0, zeroSigma );
    if ( sum( probVec ) == 0 || sum( probVec > 1e-6 ) < 5 )
        probVec = ones( size(inDistribution) ) / length( inDistribution );
    else
        probVec = probVec / sum( probVec );
    end
    probVec = cumsum( probVec );
    outDist = PerformResample( probVec, inDistribution, length( inDistribution ) );
end