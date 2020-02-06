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
    
    Tcan = params.Tcan;
    
    % In the state distribution, the first set is the internal state
    % distribution and the second is the canal distribution
    internal_x_distribution = state_distributions{1};
    internal_y_distribution = state_distributions{2};
    internal_z_distribution = state_distributions{3};
    
    canal_x_distribution = state_distributions{4};
    canal_y_distribution = state_distributions{5};
    canal_z_distribution = state_distributions{6};
    
    afferent_x_distribution = state_distributions{7};
    afferent_y_distribution = state_distributions{8};
    afferent_z_distribution = state_distributions{9};
    
    % Update canals
    act_alpha = ( Tcan * repmat( alpha, 1, length( canal_x_distribution ) ) )' + randn( length( canal_x_distribution), 3 ) * sigAlpha;
    new_canal_x_distribution = canal_x_distribution .* ( 1 - dt./tc ) + act_alpha(:,1) * dt;
    new_canal_y_distribution = canal_y_distribution .* ( 1 - dt./tc ) + act_alpha(:,2) * dt;
    new_canal_z_distribution = canal_z_distribution .* ( 1 - dt./tc ) + act_alpha(:,3) * dt;
    
    new_afferent_x_distribution = afferent_x_distribution .* ( 1 - dt./tn) - canal_x_distribution .* dt./tc + act_alpha(:,1) * dt + randn( length( afferent_z_distribution ), 1) * sigAfferent * dt ;
    new_afferent_y_distribution = afferent_y_distribution .* ( 1 - dt./tn) - canal_y_distribution .* dt./tc + act_alpha(:,2) * dt + randn( length( afferent_z_distribution ), 1) * sigAfferent * dt;
    new_afferent_z_distribution = afferent_z_distribution .* ( 1 - dt./tn) - canal_z_distribution .* dt./tc + act_alpha(:,3) * dt + randn( length( afferent_z_distribution ), 1) * sigAfferent * dt;
    
    % Resample afferents with the zero prior
    prior_afferent_x_distribution = ZeroPriorResample( new_afferent_x_distribution, sigPrior );
    %noisy_afferent_x_distribution = prior_afferent_x_distribution + randn( length( afferent_x_distribution ), 1) * sigAfferent;
    
    prior_afferent_y_distribution = ZeroPriorResample( new_afferent_y_distribution, sigPrior );
    %noisy_afferent_y_distribution = prior_afferent_y_distribution + randn( length( afferent_y_distribution ), 1) * sigAfferent;
    
    prior_afferent_z_distribution = ZeroPriorResample( new_afferent_z_distribution, sigPrior );
    %noisy_afferent_z_distribution = prior_afferent_z_distribution;
    
    % Propagate canals
%     if (tnow > 60 )
%         keyboard;
%     end
    resample_distributions = {internal_x_distribution, prior_afferent_x_distribution};
    new_internal_x_distribution = ResampleDistribution( resample_distributions, length(internal_x_distribution), 4, [10, 1] ) + sigState * randn( length( internal_x_distribution), 1 ); % * sigState * dt;
    
    resample_distributions = {internal_y_distribution, prior_afferent_y_distribution};
    new_internal_y_distribution = ResampleDistribution( resample_distributions, length(internal_y_distribution), 4, [10, 1] ) + sigState * randn( length( internal_y_distribution), 1 ); % * sigState * dt;
    
    resample_distributions = {internal_z_distribution, prior_afferent_z_distribution};
    new_internal_z_distribution = ResampleDistribution( resample_distributions, length(internal_z_distribution), 4, [10, 1] ) + sigState * randn( length( internal_z_distribution), 1 ); % * sigState * dt;
    
    new_distribution = {new_internal_x_distribution, new_internal_y_distribution, new_internal_z_distribution, ...
                        new_canal_x_distribution, new_canal_y_distribution, new_canal_z_distribution, ...
                        new_afferent_x_distribution, new_afferent_y_distribution, new_afferent_z_distribution};
    plot( tnow, mean( prior_afferent_y_distribution ), 'co', 'MarkerSize', 5, 'LineWidth', 2  );
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