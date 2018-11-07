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

function cost = FitExponentialDecay( params, t1, t2, y1, y2 )
    t1_new = t1 - t1(1);
    t2_new = t2 - t2(1);
    
    A1 = params(1);
    A2 = params(2);
    tau = params(3);
    
    ycomp1 = A1*exp( -t1_new / tau );
    ycomp2 = A2*exp( -t2_new / tau );
    
    cost = sum( ( ycomp1 - y1 ).*(ycomp1 - y1) ) + ...
           sum( ( ycomp2 - y2 ).*(ycomp2 - y2) );
end