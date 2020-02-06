%%%
% File: VelocityEstimate.m
% Author: Calvin Kuo
% Date: 10-12-2018
% Notes: This code is a restructuring of the Laurens_NatureNeuroscience_M.m
% code originally written by PAF, CJD, NK, and JSB
%
% This code provides specifically the equation for velocity estimation
% Parameters (params):
%    - gd = Gain on canal afferents
% Additional inputs (u):
%    - D = The canal afferent
%    - dirVI = Direct visual pathway
%    - VSf = Full velocity storage
% Model:
%    - VE = VSf  gd*D + dirVI

function VE = VelocityEstimate( t, u, params )
    
    % Extract parameters
    gd = params.gd;
    
    % Check to see what kind of input we have. If it is a time series, then
    % we should have a 4xn matrix input. If it is already a single angular
    % acceleration value, then we should have a 3x1 input
    
    % Case when we are passed the canal afferent at this
    % time directly
    vec_D = u.D;
    if ( size(vec_D,1) == 3 )
        assert( size(vec_D,2) == 1 );
        curr_D = vec_D;
    % Case when we are passed a time series of angular acceleration and
    % must extract the current gravito-inertial acceleration through
    % interpolation
    elseif ( size(vec_D,1) == 4 )
        tD = vec_D(1,:);
        D = vec_D(2:4,:);
        curr_D = interp1( tD, D', t )';
    % Something else, through an assertion
    else
        assert(1);
    end
    
    % Case when we are passed the full velocity storage at this time
    % directly
    vec_VSf = u.VSf;
    if ( size(vec_VSf,1) == 3 )
        assert( size(vec_VSf,2) == 1 );
        curr_VSf = vec_VSf;
    % Case when we are passed a time series of full velocity storage and
    % must extract the current full velocity storage hrough
    % interpolation
    elseif ( size(vec_VSf,1) == 4 )
        tVSf = vec_VSf(1,:);
        all_VSf = vec_VSf(2:4,:);
        
        % Interp or floor
        curr_VSf = interp1( tVSf, all_VSf', t )';
    % Something else, through an assertion
    else
        assert(1);
    end
    
    % Case when we are passed the direct visual pathway at this
    % time directly
    vec_dirVI = u.dirVI;
    if ( size(vec_dirVI,1) == 3 )
        assert( size(vec_dirVI,2) == 1 );
        curr_dirVI = vec_dirVI;
    % Case when we are passed a time series of direct visual pathway 
    % and must extract the current direct visual pathway 
    % through interpolation
    elseif ( size(vec_dirVI,1) == 4 )
        tdirVI = vec_dirVI(1,:);
        all_dirVI = vec_dirVI(2:4,:);
        
        % Interp or floor
        curr_dirVI = interp1( tdirVI, all_dirVI', t )';
    % Something else, through an assertion
    else
        assert(1);
    end
    
    % Model
    VE = curr_VSf + curr_dirVI + gd*curr_D;
end