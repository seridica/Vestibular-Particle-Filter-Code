%%%
% File: ValidateModel.m
% Author: Calvin Kuo
% Date: 10-10-2018
% Notes: This code is a restructuring of the Laurens_NatureNeuroscience_M.m
% code originally written by PAF, CJD, NK, and JSB

%% File for validating against original Laurens_Retin_vprofile_fullAdaptaion.m

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
params.tn = 100000;
params.kv = 0.247;
params.tvs = 14;

params.kf = 0.38; % Set to zero to turn off gravity
params.ts = 1/0.65;

params.go = 0; % Set to zero to turn off vision
params.ko = 0; % Set to zero to turn off vision

params.gd = 1;

params.vswitch = 0;
params.gswitch = 0;

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
mopo = 3;
if mopo == 1
    
    %%% Angular inputs from EBR paper - varies depending on plot being made
    angVMag = [0; 120*pi/180; 0];     % Velocity magnitude in rads (roll, pitch yaw)
   % Angular velocity profiles (roll, pitch, yaw) - currently uses same profile scaled to velMag
    Tmo = 316;                       % Duration of total motion profile
    t = 0:dt:Tmo;
    angPro = zeros(1,Tmo*Fs);    
    angPro(t>10&t<=110) = linspace(0,1,length(angPro(t>10&t<=110))); %angPro(t>31&t<34)=1;
    angPro(t>110&t<210)= 1;
    angPro(t>=210&t<310) = linspace(1,0,length(angPro(t>=10&t<110)));
    angVel = angVMag * angPro;
    % Angular acceleration
    angAcc = [zeros(3,1), diff(angVel,[],2)]/dt;
    t = (0:length(angAcc)-1)*dt;
    % Angular position
    angPos = cumtrapz(t,angVel,2);
    angPos(2,:) = angPos(2,:)-90*pi/180;

    %%% Translation motion profiles - acceleration
    tranAMag = [0; 0; 0];           % Translational acceleration magnitude (X,Y,Z)
    % Translational accleeration profiles
    tranPro = zeros(1,Tmo*Fs);
    tranPro(56>t&t<65) = 1;
    tranAcc = tranAMag * tranPro;
    
    % Inputs!
    u.alpha = [t;angAcc];
    u.omega = [t;angVel];
    u.acc = [t;tranAcc];
    
    % Time vector
    N = length(t);
    
elseif mopo == 2
    % Parameters for Laurens Nature paper
    % NOTE: this section was not finished.  Feel free to replace all of
    % this.
    
    kv = 0.2;                  % Cupula gain (Laurens and Angelaki)
    Tc = 4;                    % Cupula time constant (Laurens and Angelaki)
    gd = 0.8;                  % Direct vestibular gain (Laurens and Angelaki)
    Tvs = 15;                  % Velocity storage time constant (Laurens and Angelaki)

    go = 0;                    % Direct visual gain = set to zero for this
    ko = 0.0;                  % Indirect visual gain = set to zerio for this

    kf = 0.38;
    Ts = 1/0.65;
    d = 1;    % somatogravic feedback
    
    %%% Angular inputs
    % Pitch profile
    angPMag = [0; 90*pi/180; 0];
    stpT = 1.4;
    hldT = 30-stpT;
    stp = -1+2/(stpT/dt):2/(stpT/dt):1-2/(stpT/dt);
    pitPos = 10*[ones(1,10/dt) fliplr(stp) -ones(1,hldT/dt+1) stp ones(1,hldT/dt+1) fliplr(stp) -ones(1,hldT/dt+1) stp ones(1,hldT/dt+1)];
    pitVel = [0, diff(pitPos)]/dt;
    % Angular velocity profiles (roll, pitch, yaw) - currently uses same profile scaled to velMag
    angPro = [zeros(1,2/dt) (0:dt*1/0.5:1) ones(1,2/dt) fliplr(0:dt*1/0.5:1) zeros(1,10/dt)];
    angPos = angPMag * angPro;
    % Angular acceleration
    angVel = [zeros(3,1), diff(angPos,[],2)]/dt;
    angAcc = [zeros(3,1), diff(angVel,[],2)]/dt;

    %%% Translation motion profiles - acceleration
    tranAMag = [0; 0; 0];           % Translational acceleration magnitude (X,Y,Z)
    % Translational accleeration profiles
    tranPro = [zeros(1,8/dt) ones(1,2/dt) zeros(1,5.2/dt)];
    tranAcc = tranAMag * tranPro;
    
    % Time vector
    t = (0:length(angAcc)-1)*dt-10+dt;
    N = length(t);
    
    % Inputs!
    u.alpha = [t;angAcc];
    u.omega = [t;angVel];
    u.acc = [t;tranAcc];
    
elseif mopo == 3
    %%% Angular inputs from EBR paper - varies depending on plot being made
    angVMag = [0; 120*pi/180; 0];     % Velocity magnitude in rads (roll, pitch yaw)
   % Angular velocity profiles (roll, pitch, yaw) - currently uses same profile scaled to velMag
    Tmo = 200;                       % Duration of total motion profile
    t = 0:dt:Tmo;
    angPro = zeros(1,Tmo*Fs);    
    angPro(t>10&t<=11) = linspace(0,1,length(angPro(t>10&t<=11))); %angPro(t>31&t<34)=1;
    angPro(t>11&t<110)= 1;
    angPro(t>=110&t<111) = linspace(1,0,length(angPro(t>=110&t<111)));
    angVel = angVMag * angPro;
    % Angular acceleration
    angAcc = [zeros(3,1), diff(angVel,[],2)]/dt;
    t = (0:length(angAcc)-1)*dt;
    % Angular position
    angPos = cumtrapz(t,angVel,2);
    angPos(2,:) = angPos(2,:)-90*pi/180;

    %%% Translation motion profiles - acceleration
    tranAMag = [0; 0; 0];           % Translational acceleration magnitude (X,Y,Z)
    % Translational accleeration profiles
    tranPro = zeros(1,Tmo*Fs);
    tranPro(56>t&t<65) = 1;
    tranAcc = tranAMag * tranPro;
    
    % Inputs!
    u.alpha = [t;angAcc];
    u.omega = [t;angVel];
    u.acc = [t;tranAcc];
    
    % Time vector
    N = length(t);
end

% Gravity
u.grav = [0; -9.81; 0];

% Angular position (for now)
u.Tcan = [1 0 0;...
          0 1 0;...
          0 0 1];

% Nois Histories
sigC = [0.1]*pi/180;             % Canal noise (Laurens and Angelaki)
sigVS = [0.2]*pi/180;            % Vestibular storage noise (Laurens and Angelaki seemed to high, halfed it)
sigVI = [0.1]*pi/180;            % Vision noise (Laurens and Angelaki seemed high, 1/10 of listed value)

Cnoise = randn(3,length(t)) * sigC / dt;
VSnoise = randn(3,length(t)) * sigVS;
VInoise = randn(3,length(t)) * sigVI;

u.Cnoise = [t; Cnoise];
u.VSnoise = [t; VSnoise];
u.VInoise = [t; VInoise];

%% Initializing the variables
% Initializing initial conditions
C = zeros(3,1);              % Canal output signal   
D = zeros(3,1);              % Cana afferent signal
INT = zeros(3,1);            % Integral output (Endolymph)
VS = zeros(3,1);             % Velocity storage leaky integrator
VSf = zeros(3,1);            % Velocity storage full
GE = u.grav;             % Gravity estimate

%init_x = [C; D; INT; VS];
%init_store = [];
init_x = [C; D; INT; VS; VSf; GE];
%init_x = [C; D; D; INT; VS; VSf; GE];

%% Simulating the model
options = odeset('MaxStep', 0.1);
%[tinteg, full_states, store_states] = ExplicitIntegrator( @(t,y,s) SimpleLaurensVestibularModelDerivExplicit(t,y,s,u,params), t, init_x, init_store);
%[tinteg, full_states] = ode45(@(t,y) SimpleLaurensVestibularModelDeriv(t,y,u,params), t, init_x, options);
[tinteg, full_states] = ode45(@(t,y) LaurensVestibularModelDeriv(t,y,u,params), [min(t), max(t)], init_x, options);
%[tinteg, full_states] = ode45(@(t,y) CanalDynamicsDeriv(t,y,u,params), [min(t), max(t)], C);%init_x);
%[tinteg, full_states] = ode45(@(t,y) LaurensVestibularModelDeriv(t,y,u,params), [min(t), max(t)], [C;D;INT;VS;VSf]);%init_x);

%[tinteg, full_states] = ode45(@(t,y) SpineLabDualAdaptationDeriv(t,y,u,params), [min(t), max(t)], init_x);

% Get velocity estimate
nSteps = length( tinteg );
VE = zeros( 3, nSteps );
D = zeros( 3, nSteps );
for i=1:length(tinteg)
    
    % Get the fused canal signal
    u.D = full_states(i,4:6)';
    u.dirVI = [0;0;0];
    u.VSf = full_states(i,10:12)';
    VE(:,i) = VelocityEstimate( tinteg(i), u, params );
end