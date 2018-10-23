%%%
% File: RetinalSlip.m
% Author: Calvin Kuo
% Date: 10-11-2018
% Notes: This code is a restructuring of the Laurens_NatureNeuroscience_M.m
% code originally written by PAF, CJD, NK, and JSB
%
% This code provides specifically the retinal slip model. This is not an
% integrated value, but rather calculated based on other variables
% Outputs (slips):
%    - rSL = Retinal Slip
%    - indVI = indirect visual pathway
%    - dirVI = direct visual pathway
% Parameters (params):
%    - go = Direct visual gain
%    - ko = Indirect visual gain
% Additional inputs (u):
%    - angVSl = angular VSlocity (value or time signal)
%    - VS = Leaky Velocity Storage
%    - VInoise = Visual noise 
%    - Tcan = rotation matrix to canal frame
% Model:
%    - rSL = angVSl - VS + VInoise
%    - dirVI = rSL * go
%    - indVI = rSL * ko

function [rSL, indVI, dirVI] = RetinalSlip( t, u, params )
    % Extract parameters
    go = params.go;
    ko = params.ko;
    
    % Check to see what kind of input we haVS. If it is a time series, then
    % we should haVS a 4xn matrix input. If it is already a single angular
    % acceleration value, then we should haVS a 3x1 input
    
    % Case when we are passed the angular VSlocity at this time
    % directly
    VSc_omega = u.omega;
    if ( size(VSc_omega,1) == 3 )
        assert( size(VSc_omega,2) == 1 );
        curr_omega = VSc_omega;
    % Case when we are passed a time series of angular VSlocity and
    % must extract the current angular acceleration through interpolation
    elseif ( size(VSc_omega,1) == 4 )
        tomega = VSc_omega(1,:);
        omega = VSc_omega(2:4,:);
        curr_omega = interp1( tomega, omega', t )';
    % Something else, through an assertion
    else
        assert(1);
    end
    
    % Case when we are passed the VSlocity estimate at this time
    % directly
    VSc_VS = u.VS;
    if ( size(VSc_VS,1) == 3 )
        assert( size(VSc_VS,2) == 1 );
        curr_VS = VSc_VS;
    % Case when we are passed a time series of VSlocity estimate and
    % must extract the current angular acceleration through interpolation
    elseif ( size(VSc_VS,1) == 4 )
        tVS = VSc_VS(1,:);
        VS = VSc_VS(2:4,:);
        curr_VS = interp1( tVS, VS', t )';
    % Something else, through an assertion
    else
        assert(1);
    end
    
    % Case when we are passed the Tcan at this time
    % directly
    VSc_Tcan = u.Tcan;
    if ( length( size(VSc_Tcan) ) == 2 )
        Tcan = VSc_Tcan;
    % Case when we are passed a time series of Tcan and
    % must extract the current Tcan through interpolation (tricky)
    % FOR NOW - CHOOSE THE CLOSEST TCAN NEED TO FIX
    elseif ( length( size(VSc_Tcan) == 3 ) )
        Tcan_time = reshape( VSc_Tcan(1,:,:), size( VSc_Tcan, 3 ), 1 );
        tind = find( Tcan_time < t, 1, 'last' );
        Tcan = reshape( VSc_Tcan(2:4,:,tind), 3, 3 );
    % Something else, through an assertion
    else
        assert(1);
    end
    
    % Model
    rSL = AddNoise( t, Tcan*curr_omega - curr_VS, u.VInoise );
    dirVI = rSL*go;
    indVI = rSL*ko;
end