%%%
% File: KarmaliParticleDeriv.m
% Author: Calvin Kuo
% Date: 8-26-2019
% Notes: This code solves the derivative for Karmali's particle filter

function x_next = KarmaliParticleDeriv( x_pre, alpha, params, dt, st )
    
    % Extract parameters
    tc = params.tc; % Canal time constant
    t2 = params.t2;
    Q = params.Q; % input variance
    R = params.R; % afferent variance
    N = params.N; % Number of particles, for bookkeeping
    
    % Size of x_pre indicates number of particles
    assert( size( x_pre, 2 ) == N );
    canal_pre = x_pre(1:2,:); % Canals have two states
    internal_canal_pre = x_pre(3:4,:); % Internal model also has two states
    
    % Assertions for bookkeeping
    assert( length(Q) == N );
    assert( length(R) == N );
    
    % Compute variance of internal model state
    Pt = cov( internal_canal_pre' ) / dt;

    % Propagate for each particle
    internal_canal_next = zeros( size( internal_canal_pre ) );
    canal_next = zeros( size( canal_pre ) );
    vel_next = zeros( size( canal_pre ) );
    Kall = zeros( size( canal_pre ) );
    for i=1:N
        % State Dynamics
        A = [0, 1; -1/(tc(i)*t2), -(1/tc(i) + 1/t2)];
        B = [0; 1];
        C = [0, 1/t2];

        % Generate noise vectors
        input_noise = normrnd( 0, sqrt(Q(i)) );
        afferent_noise = normrnd( 0, sqrt(R(i)) );

        % Compute noisy input - Note, Karmali reports noise on the velocity
        % level, not the acceleration level
        alpha_noise = alpha; % + input_noise;

        % Canal model
        canal_deriv = A * canal_pre(:,i) + B * alpha_noise;
        canal_next(:,i) = canal_deriv * dt + canal_pre(:,i);
        in_canal = C * canal_pre(:,i) + afferent_noise;
        %Rt = std( in_canal )^2;

        % Compute "Kalman Gain"
%         if st == 1
%             K = Pt*C'/R(i);
%             %K = [0;3];
%         else
%             K = [0;3];%ones(size(R))*3; %Pt./R;
%         end

        if st == 0
            K = [0;3]; %Pt./R;
        else
            K = Pt*C'/R(i);
            %K(1) = 0;
        end
        Kall(:,i) = K;

        % Internal sensor model dynamics
        %internal_canal_next = ( K .* ( canal_next + afferent_noise ) + A .* internal_canal_pre' - vel_pre' ) ./ (K+1);
        %internal_canal_next = ( K .* ( canal_next + afferent_noise ) + A .* internal_canal_pre' - vel_pre' ) ./ (K+1);

        internal_canal_deriv = A*internal_canal_pre(:,i) + K*(in_canal - C*internal_canal_pre(:,i)) + B * input_noise;
        internal_canal_next(:,i) = internal_canal_deriv * dt + internal_canal_pre(:,i);

        % Observer
        vel_next(:,i) = K .* (in_canal - C*internal_canal_pre(:,i));
    end
        
    %x_next = [canal_next'; internal_canal_next'; vel_next'; K'];
    x_next = [canal_next; internal_canal_next; vel_next; Kall];
    %x_next = [canal_next'; internal_canal_next'; vel_next'; ones(size(R))'*(Pt./Rt)];
end