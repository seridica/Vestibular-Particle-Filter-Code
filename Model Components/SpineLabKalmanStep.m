%%%
% File: SpineLabKalmanStep.m
% Author: Calvin Kuo
% Date: 11-02-2018

function [x_t_t, sig_t_t] = SpineLabKalmanStep( x_t1_t1, sig_t1_t1, alpha, params, dt, tnow )

    % Extract parameters
    tc = params.tc;
    tn1 = params.tn1;
    tn2 = params.tn2;
    
    sigAlpha = params.sigAlpha;
    sigAfferent = params.sigAfferent;
    sigPrior = params.sigPrior;

    % State Vector - Canal State (3), Long-Term State (3), Internal Model
    % State (3)
    canal_long1_t1 = x_t1_t1(1:6);
    canal_long2_t1 = x_t1_t1(7:12);
    int_t1_t1 = x_t1_t1(13:15);
    
    % Canal update and long-term state update
    A1 = zeros(6,6);
    A1(1:3,1:3) = eye(3)*(1-dt/tc);
    A1(4:6,1:3) = eye(3)*(-dt/tc);
    A1(4:6,4:6) = eye(3)*(1-dt/tn1);
    
    B1 = [eye(3); eye(3)] * dt;
    
    C1 = [eye(3)*dt, zeros(3,3); eye(3)*dt, eye(3)];
    
    canal_long1_t = A1 * canal_long1_t1 + B1 * params.Tcan * alpha; % + C1 * [sqrt(sigAlpha) * randn(3,1); zeros(3,1)];
    
    A2 = A1;
    A2(4:6,4:6) = eye(3)*(1-dt/tn2);
    
    canal_long2_t = A2 * canal_long2_t1 + B1 * params.Tcan * alpha; % + C1 * [sqrt(sigAlpha) * randn(3,1); zeros(3,1)];
    
    % Internal state model, can update this as needed
    Aint = eye(3);
    int_t_t1 = Aint*int_t1_t1;
    sig_t_t1 = Aint'*sig_t1_t1*Aint + diag( eye(3)*abs( int_t_t1 - ( canal_long1_t(4:6) + canal_long2_t(4:6) ) / 2 ).^2 ) / 1000;
    %sig_t_t1 = diag( eye(3)*abs( int_t_t1 - ( canal_long1_t(4:6) + canal_long2_t(4:6) ) / 2 ).^2 )
    
    % Redefine a new observation distribution based on the canals and
    % previous state
    int_t_t1 = Aint*int_t1_t1; % Internal update model, currently state is just propagated through
    sig_t_t1 = Aint'*sig_t1_t1*Aint;
    
%     weights = [0.1,1,1];
    %weights = weights / sum( weights );
    canal1_est = canal_long1_t(4:6); % + sqrt(sigAfferent) * randn(3,1);
    canal2_est = canal_long2_t(4:6); % + sqrt(sigAfferent) * randn(3,1);
    
    [fs, fss, canal_weights] = CanalParticleWeighting( [canal1_est, canal2_est], [0,0,0], diag( sigPrior )' );
    weights = [100 100 100; canal_weights];
    weights(:,1) = weights(:,1) / sum(weights(:,1));
    weights(:,2) = weights(:,2) / sum(weights(:,2));
    weights(:,3) = weights(:,3) / sum(weights(:,3));
    sig_sensor = sigAlpha*dt + sigAfferent;
    combined_est = weights(1,:)'.*int_t_t1 + weights(2,:)'.*canal1_est + weights(3,:)'.*canal2_est;
    combined_sig = weights(1,:)'.*( sig_t_t1 + int_t_t1*int_t_t1' ) + ...
                   weights(2,:)'.*( sig_sensor + canal1_est*canal1_est' ) + ...
                   weights(3,:)'.*( sig_sensor + canal2_est*canal2_est' ) - combined_est*combined_est';

    keyboard;
    int_t_t = combined_est;
    sig_t_t = combined_sig;
    % Fusion
%     denom = inv( combined_sig + sigPrior );
%     
%     sig_t_t = combined_sig * sigPrior * denom;
%     int_t_t = denom * ( sigPrior*combined_est );
%             
    x_t_t = [ canal_long1_t; canal_long2_t; int_t_t ];
    
    int_back = params.Tcan'*int_t_t;
    plot( tnow, 1*int_back(2), 'ko' );
    %plot( [tnow, tnow], [int_t_t(2) - sqrt( sig_t_t(2,2) ), int_t_t(2) + sqrt( sig_t_t(2,2) )], 'k' );
    combined_back = params.Tcan'*combined_est;
    plot( tnow, 1*combined_back(2), 'bx' );
    %plot( [tnow, tnow], [combined_est(2) - sqrt( combined_sig(2,2) ), combined_est(2) + sqrt( combined_sig(2,2) )], 'b' );
    canal1_back = params.Tcan'*canal1_est;
    plot( tnow, canal1_back(2), 'rs' );
    %plot( [tnow, tnow], [canal1_est(2) - sqrt( sig_sensor(2,2) ), canal1_est(2) + sqrt( sig_sensor(2,2) )], 'r' )
    canal2_back = params.Tcan'*canal2_est;
    plot( tnow, canal2_back(2), 'rs' );
    %plot( [tnow, tnow], [canal2_est(2) - sqrt( sig_sensor(2,2) ), canal2_est(2) + sqrt( sig_sensor(2,2) )], 'r' )
    %keyboard;
end