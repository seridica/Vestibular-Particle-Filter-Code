%% PerformResample
% This is a generic helper function used by all of the resampling methods
% that performs the actual random resampling process
function outDist = PerformResample( probVec, dataVec, numSamples )
    
    % Generate random numbers for resampling
    randNums = rand( numSamples, 1 );
    if size( dataVec, 2 ) == 1
        probVec = probVec';
    end
    
    % Use random numbers to resample from the distributions
    [i,j] = find( randNums' < probVec' );
    [n,m] = unique(j);
    
    % Resampled
    tempDist = dataVec( i(m) );
    
    % Make all samples unique within the distribution
    mDist = mean( tempDist );
    stdDist = std( tempDist );
    
    % Make points unique within the resampled distribution
    uniqueDist = unique( tempDist );
    outDist = [uniqueDist; randn( length( tempDist ) - length( uniqueDist ), 1 ) * stdDist + mDist];
end