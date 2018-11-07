%%%
% File: ParticleModel.m
% Author: Calvin Kuo
% Date: 10-23-2018
% Notes: This code is a restructuring of the Laurens_NatureNeuroscience_M.m
% code originally written by PAF, CJD, NK, and JSB

%% Model parameters
% Parameters
%   - Canal Dynamics
%       : tc = Canal time constant
%       : tn = Long term adaptation time constant
%       : kv = Gain on canal afferent
%       : tvs = Velocity leakage time constant
%   - Gravity Estimation
%       : kf = Gain on gravity correction (set to 0 to turn off)
%       : ts = Gravity estimate time constant
%   - Visual Feedback
%       : go = Direct visual pathway gain (set to 0 to turn off)
%       : ko = Indirect visual pathway gain (set to 0 to turn off)
%   - Velocity Estimate
%       : gd = Gain on canal afferent

params.tc = 4;

% For dual line
params.tn1 = 40; %25; %75.9; %75;
params.tn2 = 500;

dt = 0.1;
Fs = 1/dt;

% Nois Histories
sigC = eye(3)*1.7^2; %1.7; %[3.6];%[0.1]*pi/180 / dt;             % Canal noise (Laurens and Angelaki)
sigQ = eye(3)*6.4^2; %0.01; %[14.0];

params.sigAlpha = sigC;
params.sigAfferent = sigQ;
params.sigPrior = eye(3)*1.5^2;

%% Inputs
%   - Angular position of canals (constant 3x1 or time varying 4xn)
%   - Angular velocity (constant 3x1 or time varying 4xn)
%   - Angular acceleration (constant 3x1 or time varying 4xn)
%   - Gravity (constant 3x1 or time varying 4xn)
%   - Linear acceleration (constant 3x1 or time varying 4xn)
%   - Canal Noise (constant 3x1 or time varying 4xn)
%   - Endolymph Noise Derivative (constant 3x1 or time varying 4xn)
%   - Visual Noise (constant 3x1 or time varying 4xn)
%   - Canal orientation (constant 3x3 or time varying 4x3xn)


% Motion Profile
mopo = 2;
if mopo == 1

    % Angular inputs from EBR paper - varies depending on plot being made
    angVMag = [0; 90; 0];%*pi/180; 0];     % Velocity magnitude in rads (roll, pitch yaw)
    Tmo = 220;                        % Duration of total motion profile
    t = 0:dt:Tmo;
    angPro = zeros(1,Tmo*Fs);
    angPro(t>1&t<=76) = linspace(0,1,length(angPro(t>1&t<=76)));
    angPro(t>76&t<116)= 1;
    angPro(t>=116&t<191) = linspace(1,0,length(angPro(t>=116&t<191)));
    angVel = angVMag * angPro;
    t = (0:length(angVel)-1)*dt;
    
    u.omega = [t; angVel];
    
    % Angular acceleration
    angAcc = [zeros(3,1), diff(angVel,[],2)]/dt;
    t = (0:length(angAcc)-1)*dt;
    u.alpha = [t; angAcc];
    
    % Angular position
    angPos = cumtrapz(t,angVel,2);
    angPos(2,:) = angPos(2,:)-90; %*pi/180;

    % Translation motion profiles - acceleration
    tranAMag = [0; 0; 0];
    % Translational accleeration profiles
    tranPro = zeros(1,Tmo*Fs);
    tranPro(56>t&t<65) = 1;
    tranAcc = tranAMag * tranPro;
    
    u.acc = [t; tranAcc];
    
elseif mopo == 2
    % Parameters for Laurens Nature paper
    % NOTE: this section was not finished.  Feel free to replace all of
    % this.
    
    % Angular inputs from EBR paper - varies depending on plot being made
    angVMag = [0; 90; 0]; %*pi/180; 0];     % Velocity magnitude in rads (roll, pitch yaw)
    Tmo = 80                        % Duration of total motion profile
    t = 0:dt:Tmo;
    angPro = zeros(1,Tmo*Fs);
    angPro(t>0.5&t<=1.5) = linspace(0,1,length(angPro(t>0.5&t<=1.5)));
    angPro(t>1.5&t<46.5)= 1;
    angPro(t>=46.5&t<47.5) = linspace(1,0,length(angPro(t>=46.5&t<47.5)));
    %angPro(t>1.5&t<101.5) = 1;
    %angPro(t>=101.5&t<102.5) = linspace(1,0,length(angPro(t>=101.5&t<102.5)));
    angVel = angVMag * angPro;
    t = (0:length(angVel)-1)*dt;
    
    u.omega = [t; angVel];
    
    % Angular acceleration
    angAcc = [zeros(3,1), diff(angVel,[],2)]/dt;
    t = (0:length(angAcc)-1)*dt;
    u.alpha = [t; angAcc];
    
    % Angular position
    angPos = cumtrapz(t,angVel,2);
    angPos(2,:) = angPos(2,:)-90; %*pi/180;

    % Translation motion profiles - acceleration
    tranAMag = [0; 0; 0];
    % Translational accleeration profiles
    tranPro = zeros(1,Tmo*Fs);
    tranPro(56>t&t<65) = 1;
    tranAcc = tranAMag * tranPro;
    
    u.acc = [t; tranAcc];

end

% Gravity
u.grav = [0; -9.81; 0];

% Angular position (for now)
params.Tcan = [1 0 0;...
          0 1 0;...
          0 0 1];
      
% params.Tcan = [1 0 0; ...
%                0 cos(pi/4) -sin(pi/4); ...
%                0 sin(pi/4) cos(pi/4)];

% params.Tcan = [1 0 0; ...
%                0 cos(pi/4) -sin(pi/4); ...
%                0 sin(pi/4) cos(pi/4)] * ...
%               [cos(pi/4) -sin(pi/4) 0; ...
%                sin(pi/4) cos(pi/4) 0; ...
%                0 0 1];

%% Initializing the variables
% Initializing initial conditions
C = zeros(3,1);              % Canal output signal   
D = zeros(3,1);              % Cana afferent signal
INT = zeros(3,1);            % Integral output (Endolymph)
VS = zeros(3,1);             % Velocity storage leaky integrator
VSf = zeros(3,1);            % Velocity storage full
GE = u.grav;             % Gravity estimate

%init_x = [C; D; INT; VS; VSf; GE];
%% Simulating with Spine Lab Kalman Filter
figure(1); clf; hold on;
x_1_1 = zeros(15,1);
sig_1_1 = params.sigAfferent;

states = zeros( 15, length(t) );
sigs = zeros( 3, length(t) );

states(:,1) = x_1_1;
sigs(:,1) = diag( sig_1_1 );

x_prev = x_1_1;
sig_prev = sig_1_1;

plot( t(1), x_1_1(2), 'ko' );
plot( [t(1), t(1)], [x_1_1(2) - sig_1_1(2), x_1_1(2) + sig_1_1(2)], 'k' );

for i=2:length( t )
    alpha = angAcc(:,i);
    [x_new, sig_new] = SpineLabKalmanStep( x_prev, sig_prev, alpha, params, dt, t(i) );
    
    states(:,i) = x_new;
    sigs(:,i) = diag( sig_new );
    
    x_prev = x_new;
    sig_prev = sig_new;
end

if mopo == 2
    % Positive and Negative Curve Fit
    full_rot = params.Tcan'*states(13:15,:);
    [m,i] = max( full_rot(2,:) );
    t1 = t(i:find( t > 19, 1, 'first' ) );
    y1 = full_rot(2,i:find( t > 19, 1, 'first' ) );
    
    [m,i] = min( full_rot(2,:) );
    t2 = t(i:end );
    y2 = full_rot(2,i:end );
    
    x0 = [y1(1), y2(1), 30];
    x = fminsearch(@(x) FitExponentialDecay( x, t1, t2, y1, y2 ), x0 )
end