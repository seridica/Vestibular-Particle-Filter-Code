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
params.tn_internal = 75;
params.kv = 0.52; %0.247;
params.tvs = 7.7; %14;

% Switches
params.vswitch = 0;
params.gswitch = 0;

% For dual line
params.tn1 = 40; %25; %75.9; %75;
params.tn2 = 10000;
params.sigprior = ones(1,3)*6^2; %0.005;

params.kf = 0.38; % Set to zero to turn off gravity
params.ts = 1/0.65;

params.go = 0; % Set to zero to turn off vision
params.ko = 0; % Set to zero to turn off vision

params.gd = 1;

dt = 0.1;
Fs = 1/dt;

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
    Tmo = 200; %200; %80                        % Duration of total motion profile
    t = 0:dt:Tmo;
    angPro = zeros(1,Tmo*Fs);
    angPro(t>0.5&t<=1.5) = linspace(0,1,length(angPro(t>0.5&t<=1.5)));
    %angPro(t>1.5&t<46.5)= 1;
    %angPro(t>=46.5&t<47.5) = linspace(1,0,length(angPro(t>=46.5&t<47.5)));
    angPro(t>1.5&t<101.5) = 1;
    angPro(t>=101.5&t<102.5) = linspace(1,0,length(angPro(t>=101.5&t<102.5)));
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
u.Tcan = [1 0 0;...
          0 1 0;...
          0 0 1];

% Nois Histories
sigC = [0.1]; %[1.7]; %[3.6];%[0.1]*pi/180 / dt;             % Canal noise (Laurens and Angelaki)
sigQ = [0.1];%[6.4]; %[14.0];

sigVS = [0.2]*pi/180;            % Vestibular storage noise (Laurens and Angelaki seemed to high, halfed it)
sigVI = [0.1]*pi/180;            % Vision noise (Laurens and Angelaki seemed high, 1/10 of listed value)

C1noise = randn(3,length(t)) * sigC; % / dt;
C2noise = randn(3,length(t)) * sigC; % / dt;
VSnoise = randn(3,length(t)) * sigVS / dt;
VInoise = randn(3,length(t)) * sigVI / dt;
dVSnoise = [zeros(3,1), diff(VSnoise,[],2)];

u.Qnoise = eye(3)*sigQ;
u.Cnoise = eye(3)*sigC; %[t; C1noise];
u.dVSnoise = [t; dVSnoise];
u.VInoise = [t; VInoise];

%% Initializing the variables
% Initializing initial conditions
C = zeros(3,1);              % Canal output signal   
D = zeros(3,1);              % Cana afferent signal
INT = zeros(3,1);            % Integral output (Endolymph)
VS = zeros(3,1);             % Velocity storage leaky integrator
VSf = zeros(3,1);            % Velocity storage full
GE = u.grav;             % Gravity estimate

%init_x = [C; D; INT; VS; VSf; GE];
%% Simulating with Karmali Derivative
nChannels = 10;
init_x = zeros(nChannels*3*4, 1);%[randn(nChannels*2*3,1)*sigVI];
init_store = [0;0;0;sigQ;sigQ;sigQ;0;0;0;randn(nChannels*6*3,1)*sigC];
options = odeset('MaxStep', 0.1);
[tinteg, full_states, full_store] = ExplicitIntegrator(@(t,y,s) KarmaliDerivExplicit(t,y,s,u,params), t, init_x, init_store);

%% Simulating with my velocity storage model
% nChannels = 10;
% init_x = zeros( nChannels*3*8, 1 );
% init_store = [0;0;0;sigQ;sigQ;sigQ;randn(nChannels*2*3,1)];
% [integ, full_states, full_store] = ExplicitIntegrator(@(t,y,s) VelocityStorageFilterDerivExplicit( t,y,s,u,params ), t, init_x, init_store );