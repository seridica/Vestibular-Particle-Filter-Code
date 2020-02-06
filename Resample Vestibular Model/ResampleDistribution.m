%%%
% File: ResampleDistribution.m
% Author: Calvin Kuo
% Date: 02-22-2019
%
% Notes - This code performs the resampling operation and returns the new
% sample distribution. There are several different types of resampling
% procedures that can be performed here.

function outDist = ResampleDistribution(inDists, numSamples, funcType, funcParams)

    % Use the default resample method
    if nargin == 2
        funcType = 1;
        funcParams = [];
    % Use the default parameters for each resample method
    elseif nargin == 3
        switch funcType
            case 1
                funcParams = [];
            case 2
                nDists = length( inDists );
                funcParams = [ones( 1, nDists ) / nDists];
            case 3
                funcParams = [];
            case 4
                nDists = length( inDists );
                funcParams = [ones( 1, nDists ) / nDists];
        end
    end
    
    % Resampler methods
    switch funcType
        case 1
            outDist = UniformResampling( inDists, numSamples, funcParams );
        case 2
            outDist = WeightedResampling( inDists, numSamples, funcParams );
        case 3
            outDist = BayesianOptimalResampling( inDists, numSamples, funcParams );
        case 4
            outDist = BayesianWeightedResampling( inDists, numSamples, funcParams );
    end
end

%% UniformResampling
% This method uniformally resamples the input distributions.
function outDist = UniformResampling( inDists, numSamples, params )
    numDists = length( inDists );
    
    % FIgure out total number of samples
    totalSamples = 0;
    for i=1:numDists
        thisDist = inDists{i};
        totalSamples = totalSamples + length( thisDist );
    end
    
    % Set up probability and data vectors
    dataVec = zeros( totalSamples, 1 );
    currInd = 1;
    for i=1:numDists
        thisDist = inDists{i};
        nSamples = length( thisDist );
        
        if size( thisDist, 1 ) == 1
            dataVec(currInd:(currInd + nSamples - 1)) = thisDist';
        else
            dataVec(currInd:(currInd + nSamples - 1)) = thisDist;
        end
        
        currInd = currInd + nSamples;
    end
    
    probVec = ones( totalSamples, 1 ) / totalSamples;
    
    assert( length( probVec ) == length( dataVec ) );
    
    % Add noise to the output distribution
    outDist = PerformResample( probVec, dataVec, numSamples );
end

%% WeightedResampling
% This method applies fixed weights to the input distributions (can mimic
% uniform resampling if all weights are equal).
%
% Params here are RELATIVE weights of the n distributions
function outDist = WeightedResampling( inDists, numSamples, params )
    assert( length( params ) == length( inDists ) );
    numDists = length( inDists );
    
    % Normalize distribution relative weights
    relWeights = params / sum( params );
    
    % Figure out total number of samples
    totalSamples = 0;
    for i=1:numDists
        thisDist = inDists{i};
        totalSamples = totalSamples + length( thisDist );
    end
    
    % Calculate weights
    dataVec = zeros( totalSamples, 1 );
    probVec = ones( totalSamples, 1 ) / totalSamples;
    currInd = 1;
    for i=1:numDists
        thisDist = inDists{i};
        nSamples = length( thisDist );
        
        if size( thisDist, 1 ) == 1
            dataVec(currInd:(currInd + nSamples - 1)) = thisDist';
        else
            dataVec(currInd:(currInd + nSamples - 1)) = thisDist;
        end
        
        probVec(currInd:(currInd + nSamples - 1)) = probVec(currInd:(currInd + nSamples - 1)) * relWeights(i);
        
        currInd = currInd + nSamples;
    end
    
    % Check
    probVec = probVec / sum( probVec );
    probVec = cumsum( probVec );
    outDist = PerformResample( probVec, dataVec, numSamples );
end

%% BayesianOptimalResampling
% This method applies weights on the distributions according to their
% vairances (as in a Bayesian optimal approach)
function outDist = BayesianOptimalResampling( inDists, numSamples, params )
    numDists = length( inDists );
    
    % Initial weights
    distVars = zeros( 1, numDists );
    
    % Figure out total number of samples, and also relative Bayesian
    % weights
    totalSamples = 0;
    for i=1:numDists
        thisDist = inDists{i};
        totalSamples = totalSamples + length( thisDist );
        
        distVars(i) = std( thisDist )^2;
    end
    
    %% TODO - DO THIS BETTER
    if numDists == 2
        relWeights = [distVars(2) / sum( distVars ), distVars(1) / sum( distVars )];
    elseif numDists == 3
        denom = distVars(1)*distVars(2) + distVars(2)*distVars(3) + distVars(1)*distVars(3);
        relWeights = [distVars(2)*distVars(3), distVars(1)*distVars(3), distVars(1)*distVars(2)] / denom;
    else
        assert(0);
    end
    
    % Calculate weights
    dataVec = zeros( totalSamples, 1 );
    probVec = ones( totalSamples, 1 ) / totalSamples;
    currInd = 1;
    for i=1:numDists
        thisDist = inDists{i};
        nSamples = length( thisDist );
        
        if size( thisDist, 1 ) == 1
            dataVec(currInd:(currInd + nSamples - 1)) = thisDist';
        else
            dataVec(currInd:(currInd + nSamples - 1)) = thisDist;
        end
        
        probVec(currInd:(currInd + nSamples - 1)) = probVec(currInd:(currInd + nSamples - 1)) * relWeights(i);
        
        currInd = currInd + nSamples;
    end
    
    % Check
    %relWeights
    %sum( probVec )
    %assert( abs( sum( probVec ) - 1 ) < 1e-6 );
    probVec = probVec / sum( probVec );
    probVec = cumsum( probVec );
    outDist = PerformResample( probVec, dataVec, numSamples );
end

%% BayesianWeightedResampling
% This is a modification on the Bayesian Optimal Resampling with an
% additional bias on the Bayesian weights.
function outDist = BayesianWeightedResampling( inDists, numSamples, params )
    numDists = length( inDists );
    
    % Initial weights
    distVars = zeros( numDists, 1 );
    
    % Normalize distribution relative weights
    relBiases = params / sum( params );
    
    % Figure out total number of samples
    totalSamples = 0;
    for i=1:numDists
        thisDist = inDists{i};
        totalSamples = totalSamples + length( thisDist );
        
        distVars(i) = std( thisDist )^2;
    end
    
    %% TODO - DO THIS BETTER
    if numDists == 2
        relWeights = [distVars(2) / sum( distVars ), distVars(1) / sum( distVars )];
    elseif numDists == 3
        denom = distVars(1)*distVars(2) + distVars(2)*distVars(3) + distVars(1)*distVars(3);
        relWeights = [distVars(2)*distVars(3), distVars(1)*distVars(3), distVars(1)*distVars(2)] / denom;
    else
        assert(0);
    end
    
    % Calculate weights
    dataVec = zeros( totalSamples, 1 );
    probVec = ones( totalSamples, 1 ) / totalSamples;
    currInd = 1;
    for i=1:numDists
        thisDist = inDists{i};
        nSamples = length( thisDist );
        
        if size( thisDist, 1 ) == 1
            dataVec(currInd:(currInd + nSamples - 1)) = thisDist';
        else
            dataVec(currInd:(currInd + nSamples - 1)) = thisDist;
        end
        
        probVec(currInd:(currInd + nSamples - 1)) = probVec(currInd:(currInd + nSamples - 1)) * relBiases(i) * relWeights(i);
        
        currInd = currInd + nSamples;
    end
    
    % Check
    %relWeights
    %sum( probVec )
    %assert( abs( sum( probVec ) - 1 ) < 1e-6 );
    probVec = probVec / sum( probVec );
    probVec = cumsum( probVec );
    outDist = PerformResample( probVec, dataVec, numSamples );
end

