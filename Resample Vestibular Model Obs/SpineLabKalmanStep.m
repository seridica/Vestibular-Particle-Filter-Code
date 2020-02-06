%%%
% File: SpineLabKalmanStep.m
% Author: Calvin Kuo
% Date: 11-02-2018

function [x_t_t, sig_t_t] = SpineLabKalmanStep( x_t1_t1, sig_t1_t1, alpha, params, dt, tnow )
    [x_t_t, sig_t_t] = SpineLabKalmanStep1( x_t1_t1, sig_t1_t1, alpha, params, dt, tnow );
end

function [x_t_t, sig_t_t] = SpineLabKalmanStep1( x_t1_t1, sig_t1_t1, alpha, params, dt, tnow )

    % Extract parameters
    tc1 = params.tc1;
    tc2 = params.tc2;
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
    A1(1:3,1:3) = eye(3)*(1-dt/tc1);
    A1(4:6,1:3) = eye(3)*(-dt/tc1);
    A1(4:6,4:6) = eye(3)*(1-dt/tn1);
    
    B1 = [eye(3); eye(3)] * dt;
    
    C1 = [eye(3)*dt, zeros(3,3); eye(3)*dt, eye(3)];
    
    canal_long1_t = A1 * canal_long1_t1 + B1 * params.Tcan * alpha; % + C1 * [sqrt(sigAlpha) * randn(3,1); zeros(3,1)];
    
    A2 = zeros(6,6);
    A2(1:3,1:3) = eye(3)*(1-dt/tc2);
    A2(4:6,1:3) = eye(3)*(-dt/tc2);
    A2(4:6,4:6) = eye(3)*(1-dt/tn2);
    
    canal_long2_t = A2 * canal_long2_t1 + B1 * params.Tcan * alpha; % + C1 * [sqrt(sigAlpha) * randn(3,1); zeros(3,1)];
    
    canal1_est = canal_long1_t(4:6); % + sqrt(sigAfferent) * randn(3,1);
    canal2_est = canal_long2_t(4:6); % + sqrt(sigAfferent) * randn(3,1);
    
    sig_sensor = sigAlpha*dt*dt + sigAfferent;
    [fs, fss, canal_weights] = CanalParticleWeighting( [canal1_est, canal2_est], [0,0,0], [diag( sigPrior )'] );
    new_canal_est = fs;
    new_canal_sig = sig_sensor + diag( canal_weights(1,:)'.*(canal1_est.*canal1_est) ) + diag( canal_weights(2,:)'.*( canal2_est.*canal2_est ) ) - diag( new_canal_est.*new_canal_est );
    
    % Internal state model, can update this as needed
    Aint = eye(3);
    int_t_t1 = Aint*int_t1_t1; % Internal update model, currently state is just propagated through
    sig_t_t1 = ( Aint'*sig_t1_t1*Aint ); %+ sigAfferent*dt );
    
    fss_mat = new_canal_sig / dt;
    %weights = [[10, 10, 10]*sig_sensor*inv(sig_sensor+sig_t_t1); canal_weights*sig_t_t1*inv(sig_sensor+sig_t_t1)];
    weights = [[100,100,100]*fss_mat; [1,1,1]*sig_t_t1];
    weights(:,1) = weights(:,1) / sum(weights(:,1));
    weights(:,2) = weights(:,2) / sum(weights(:,2));
    weights(:,3) = weights(:,3) / sum(weights(:,3));
    combined_est = weights(1,:)'.*int_t_t1 + weights(2,:)'.*fs;
    combined_sig = weights(1,:)'.*( sig_t_t1 + diag(int_t_t1.*int_t_t1) ) + ...
                   weights(2,:)'.*( fss_mat + diag(fs.*fs) ) - diag(combined_est.*combined_est);
    
    int_t_t = combined_est;
    sig_t_t = combined_sig; %[fss_mat(2,2), 0, 0; 0, combined_sig(2,2), 0; 0, 0, 1];%;
    
    % Fusion
%     denom = inv( combined_sig + sigPrior );
%     
%     sig_t_t = combined_sig * sigPrior * denom;
%     int_t_t = denom * ( sigPrior*combined_est );
%             
    x_t_t = [ canal_long1_t; canal_long2_t; int_t_t ];
    %keyboard;
    int_back = params.Tcan'*int_t_t;
    sig_back = params.Tcan'*sig_t_t*params.Tcan;
    plot( tnow, int_back(2), 's', 'Color', [0.0, 0.0, 0.7] );
    plot( [tnow, tnow], [int_back(2) - sqrt( sig_back(2,2) ), int_back(2) + sqrt( sig_back(2,2) )], 'Color', [0.3, 0.3, 1.0] );
    %combined_back = params.Tcan'*combined_est;
    %plot( tnow, 1*combined_back(2), 'bx' );
    %plot( [tnow, tnow], [combined_est(2) - sqrt( combined_sig(2,2) ), combined_est(2) + sqrt( combined_sig(2,2) )], 'b' );
%     canal1_back = params.Tcan'*canal1_est;
     sensor_back = params.Tcan'*sig_sensor*params.Tcan;
     canal_back = params.Tcan'*new_canal_sig*params.Tcan;
%     plot( tnow, canal1_back(2), 's', 'Color', [0.7, 0.7, 0.7],'HandleVisibility','off'  );
%     plot( [tnow, tnow], [canal1_back(2) - sqrt( sensor_back(2,2) ), canal1_back(2) + sqrt( sensor_back(2,2) )], 'r' )
%     canal2_back = params.Tcan'*canal2_est;
%     plot( tnow, canal2_back(2), 's', 'Color', [0.7, 0.7, 0.7],'HandleVisibility','off'  );
%     plot( [tnow, tnow], [canal2_back(2) - sqrt( sensor_back(2,2) ), canal2_back(2) + sqrt( sensor_back(2,2) )], 'r' )
    fr = params.Tcan' * fs;
    plot( tnow, fr(2), 'kx' );
    plot( [tnow, tnow], [fr(2) - sqrt( canal_back(2,2) ), fr(2) + sqrt( canal_back(2,2) )], 'Color', [0.7, 0.7, 0.7] )
end

function [x_t_t, sig_t_t] = SpineLabKalmanStep2( x_t1_t1, sig_t1_t1, alpha, params, dt, tnow )

    % Extract parameters
    tc = params.tc;
    tn1 = params.tn1;
    tn2 = params.tn2;
    
    sigAlpha = params.sigAlpha;
    sigAfferent = params.sigAfferent;
    sigPrior = params.sigPrior;
    sigState = params.sigState;

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
    
    canal1_est = canal_long1_t(4:6); % + sqrt(sigAfferent) * randn(3,1);
    canal2_est = canal_long2_t(4:6); % + sqrt(sigAfferent) * randn(3,1);
    
    sig_sensor = ( sigAlpha + sigAfferent ) * dt;
    [fs, fss, canal_weights] = CanalParticleWeighting( [canal1_est, canal2_est], [0,0,0], [diag( sigPrior )'] );
    new_canal_est = fs;
    new_canal_sig = sig_sensor + diag( canal_weights(1,:)'.*(canal1_est.*canal1_est) ) + diag( canal_weights(2,:)'.*( canal2_est.*canal2_est ) ) - diag( new_canal_est.*new_canal_est );
    
    % Internal state model, can update this as needed
    Aint = eye(3);
    int_t_t1 = Aint*int_t1_t1; % Internal update model, currently state is just propagated through
    sig_t_t1 = ( Aint'*sig_t1_t1*Aint + sigState * dt );
    
    fss_mat = new_canal_sig;
    %weights = [[10, 10, 10]*sig_sensor*inv(sig_sensor+sig_t_t1); canal_weights*sig_t_t1*inv(sig_sensor+sig_t_t1)];
    weights = [[(5000),(5000),(5000)]*fss_mat; [1,1,1]*sig_t_t1];
    weights(:,1) = weights(:,1) / sum(weights(:,1));
    weights(:,2) = weights(:,2) / sum(weights(:,2));
    weights(:,3) = weights(:,3) / sum(weights(:,3));
    combined_est = weights(1,:)'.*int_t_t1 + weights(2,:)'.*fs;
    combined_sig = weights(1,:)'.*( sig_t_t1 + diag(int_t_t1.*int_t_t1) ) + ...
                   weights(2,:)'.*( fss_mat + diag(fs.*fs) ) - diag(combined_est.*combined_est);
    
    int_t_t = combined_est;
    sig_t_t = combined_sig; %[fss_mat(2,2), 0, 0; 0, combined_sig(2,2), 0; 0, 0, 1];%;
    
    % Fusion
%     denom = inv( combined_sig + sigPrior );
%     
%     sig_t_t = combined_sig * sigPrior * denom;
%     int_t_t = denom * ( sigPrior*combined_est );
%             
    x_t_t = [ canal_long1_t; canal_long2_t; int_t_t ];
    %keyboard;
    int_back = params.Tcan'*int_t_t;
    sig_back = params.Tcan'*sig_t_t*params.Tcan;
    plot( tnow, int_back(2), 's', 'Color', [0.0, 0.0, 0.7] );
    plot( [tnow, tnow], [int_back(2) - sqrt( sig_back(2,2) ), int_back(2) + sqrt( sig_back(2,2) )], 'Color', [0.3, 0.3, 1.0] );
    %combined_back = params.Tcan'*combined_est;
    %plot( tnow, 1*combined_back(2), 'bx' );
    %plot( [tnow, tnow], [combined_est(2) - sqrt( combined_sig(2,2) ), combined_est(2) + sqrt( combined_sig(2,2) )], 'b' );
%     canal1_back = params.Tcan'*canal1_est;
     sensor_back = params.Tcan'*sig_sensor*params.Tcan;
     canal_back = params.Tcan'*new_canal_sig*params.Tcan;
%     plot( tnow, canal1_back(2), 's', 'Color', [0.7, 0.7, 0.7],'HandleVisibility','off'  );
%     plot( [tnow, tnow], [canal1_back(2) - sqrt( sensor_back(2,2) ), canal1_back(2) + sqrt( sensor_back(2,2) )], 'r' )
%     canal2_back = params.Tcan'*canal2_est;
%     plot( tnow, canal2_back(2), 's', 'Color', [0.7, 0.7, 0.7],'HandleVisibility','off'  );
%     plot( [tnow, tnow], [canal2_back(2) - sqrt( sensor_back(2,2) ), canal2_back(2) + sqrt( sensor_back(2,2) )], 'r' )
    fr = params.Tcan' * fs;
    plot( tnow, fr(2), 'kx' );
    plot( [tnow, tnow], [fr(2) - sqrt( canal_back(2,2) ), fr(2) + sqrt( canal_back(2,2) )], 'Color', [0.7, 0.7, 0.7] )
end



function [x_t_t, sig_t_t] = SpineLabKalmanStep3( x_t1_t1, sig_t1_t1, alpha, params, dt, tnow )

    % Extract parameters
    tc1 = params.tc1;
    tc2 = params.tc2;
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
    A1(1:3,1:3) = eye(3)*(1-dt/tc1);
    A1(4:6,1:3) = eye(3)*(-dt/tc1);
    A1(4:6,4:6) = eye(3)*(1-dt/tn1);
    
    B1 = [eye(3); eye(3)] * dt;
    
    C1 = [eye(3)*dt, zeros(3,3); eye(3)*dt, eye(3)];
    
    canal_long1_t = A1 * canal_long1_t1 + B1 * params.Tcan * alpha; % + C1 * [sqrt(sigAlpha) * randn(3,1); zeros(3,1)];
    
    A2 = zeros(6,6);
    A2(1:3,1:3) = eye(3)*(1-dt/tc2);
    A2(4:6,1:3) = eye(3)*(-dt/tc2);
    A2(4:6,4:6) = eye(3)*(1-dt/tn2);
    
    canal_long2_t = A2 * canal_long2_t1 + B1 * params.Tcan * alpha; % + C1 * [sqrt(sigAlpha) * randn(3,1); zeros(3,1)];
    
    canal1_est = canal_long1_t(4:6); % + sqrt(sigAfferent) * randn(3,1);
    canal2_est = canal_long2_t(4:6); % + sqrt(sigAfferent) * randn(3,1);
    
    sig_sensor = sigAlpha*dt*dt + sigAfferent / dt;
    [fs, fss, canal_weights] = CanalParticleWeighting( [canal1_est, canal2_est], [0,0,0], [diag( sigPrior )'] );
    new_canal_est = fs;
    new_canal_sig = sig_sensor + diag( canal_weights(1,:)'.*(canal1_est.*canal1_est) ) + diag( canal_weights(2,:)'.*( canal2_est.*canal2_est ) ) - diag( new_canal_est.*new_canal_est );
    
    % Internal state model, can update this as needed
    Aint = eye(3);
    %int_t_t1 = Aint*int_t1_t1; % Internal update model, currently state is just propagated through
    %sig_t_t1 = ;
    
    %fss_mat = new_canal_sig;
    %weights = [[10, 10, 10]*sig_sensor*inv(sig_sensor+sig_t_t1); canal_weights*sig_t_t1*inv(sig_sensor+sig_t_t1)];
    %weights = [[100,100,100]*fss_mat; [1,1,1]*sig_t_t1];
    %weights(:,1) = weights(:,1) / sum(weights(:,1));
    %weights(:,2) = weights(:,2) / sum(weights(:,2));
    %weights(:,3) = weights(:,3) / sum(weights(:,3));
    %combined_est = weights(1,:)'.*int_t_t1 + weights(2,:)'.*fs;
    %combined_sig = weights(1,:)'.*( sig_t_t1 + diag(int_t_t1.*int_t_t1) ) + ...
    %               weights(2,:)'.*( fss_mat + diag(fs.*fs) ) - diag(combined_est.*combined_est);
    
    %int_t_t = combined_est;
    %sig_t_t = combined_sig; %[fss_mat(2,2), 0, 0; 0, combined_sig(2,2), 0; 0, 0, 1];%;
    %sig_t1_t = sig_t1_t1 + sigAfferent*dt;
    
    K1gain = sig_t1_t1*inv( sig_t1_t1 + new_canal_sig );
    %K2gain = 1000/dt*new_canal_sig*inv( sig_t1_t + 1000/dt*new_canal_sig );
    int_t_t = int_t1_t1 + K1gain * ( fs - int_t1_t1 );
    %sig_t_t = ( eye(3) - K1gain ) * sig_t1_t1 + sigAfferent*dt; 
    sig_t_t = ( eye(3) - K1gain ) * sig_t1_t1 + dt * sigAfferent; %* ( diag(fs.*fs) + diag( int_t1_t1.*int_t1_t1 ) - diag( int_t_t.*int_t_t ) );
    %sigAfferent * dt;
    %dt*dt*( diag(fs.*fs) + diag( int_t1_t1 .* int_t1_t1 ) - diag( int_t_t .* int_t_t ) );%sigAfferent*dt;
    %sig_t_t = K1gain * ( new_canal_sig + diag(fs.*fs) ) + ...
    %        K2gain * ( sig_t1_t + diag( int_t1_t1 .* int_t1_t1 ) ) - ...
    %        diag( int_t_t .* int_t_t );
    
    %int_t_t = int_t1_t1 + ( sig_t1_t1 * inv( new_canal_sig ) ) * dt * ( fs - int_t1_t1 );
    %sig_t_t = sig_t1_t1 * (1+2*dt) - sig_t1_t1 * inv( new_canal_sig ) * sig_t1_t1 * dt + sigAfferent * dt;
    
    %sig_t_t = (new_canal_sig - sig_t1_t1*dt)*inv(new_canal_sig)*( sig_t1_t1 + diag(int_t1_t1.*int_t1_t1) ) + ...
    %          (sig_t1_t1*dt)*inv(new_canal_sig)*( new_canal_sig + diag(fs.*fs) ) - diag(int_t_t.*int_t_t);
    %if (tnow > 0.5)
    %   keyboard;
    %end
    
    x_t_t = [ canal_long1_t; canal_long2_t; int_t_t ];
    %keyboard;
    int_back = params.Tcan'*int_t_t;
    sig_back = params.Tcan'*sig_t_t*params.Tcan;
    plot( tnow, int_back(2), 'gs' );
    plot( [tnow, tnow], [int_back(2) - sqrt( sig_back(2,2) ), int_back(2) + sqrt( sig_back(2,2) )], 'k' );
    %combined_back = params.Tcan'*combined_est;
    %plot( tnow, 1*combined_back(2), 'bx' );
    %plot( [tnow, tnow], [combined_est(2) - sqrt( combined_sig(2,2) ), combined_est(2) + sqrt( combined_sig(2,2) )], 'b' );
    canal1_back = params.Tcan'*canal1_est;
    sensor_back = params.Tcan'*sig_sensor*params.Tcan;
    plot( tnow, canal1_back(2), 's', 'Color', [0.7, 0.7, 0.7],'HandleVisibility','off'  );
    plot( [tnow, tnow], [canal1_back(2) - sqrt( sensor_back(2,2) ), canal1_back(2) + sqrt( sensor_back(2,2) )], 'r' )
    canal2_back = params.Tcan'*canal2_est;
    plot( tnow, canal2_back(2), 's', 'Color', [0.7, 0.7, 0.7],'HandleVisibility','off'  );
    plot( [tnow, tnow], [canal2_back(2) - sqrt( sensor_back(2,2) ), canal2_back(2) + sqrt( sensor_back(2,2) )], 'r' )
    fr = params.Tcan' * fs;
    plot( tnow, fr(2), 'bx' );
end