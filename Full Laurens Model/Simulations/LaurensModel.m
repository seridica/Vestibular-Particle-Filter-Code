%%%
% File: LaurensModel.m
% Author: Calvin Kuo
% Date: 10-10-2018
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

params.tc = 7.3; %4;
params.tn = 25; %75;
params.kv = 0.247; %0.52; %0.247;
params.tvs = 7.7; %14;

% Switches
params.vswitch = 0;
params.gswitch = 0;

% For dual line
params.tn1 = 10000; %75.9; %75;
params.tn2 = 10000;
params.sigprior = ones(1,3)*0.5;

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
mopo = 1;
if mopo == 1

    % Angular inputs from EBR paper - varies depending on plot being made
    angVMag = [0; 90*pi/180; 0];     % Velocity magnitude in rads (roll, pitch yaw)
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
    angPos(2,:) = angPos(2,:)-90*pi/180;

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
    angVMag = [0; 90*pi/180; 0];     % Velocity magnitude in rads (roll, pitch yaw)
    Tmo = 220;                        % Duration of total motion profile
    t = 0:dt:Tmo;
    angPro = zeros(1,Tmo*Fs);
    angPro(t>1&t<=2) = linspace(0,1,length(angPro(t>1&t<=2)));
    angPro(t>2&t<100)= 1;
    angPro(t>=100&t<101) = linspace(1,0,length(angPro(t>=100&t<101)));
    angVel = angVMag * angPro;
    t = (0:length(angVel)-1)*dt;
    
    u.omega = [t; angVel];
    
    % Angular acceleration
    angAcc = [zeros(3,1), diff(angVel,[],2)]/dt;
    t = (0:length(angAcc)-1)*dt;
    u.alpha = [t; angAcc];
    
    % Angular position
    angPos = cumtrapz(t,angVel,2);
    angPos(2,:) = angPos(2,:)-90*pi/180;

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
sigC = [0.1]*pi/180;             % Canal noise (Laurens and Angelaki)
sigVS = [0.2]*pi/180;            % Vestibular storage noise (Laurens and Angelaki seemed to high, halfed it)
sigVI = [0.1]*pi/180;            % Vision noise (Laurens and Angelaki seemed high, 1/10 of listed value)

Cnoise = randn(3,length(t)) * sigC / dt;
VSnoise = randn(3,length(t)) * sigVS / dt;
VInoise = randn(3,length(t)) * sigVI / dt;
dVSnoise = [zeros(3,1), diff(VSnoise,[],2)];

u.Cnoise = [t; Cnoise];
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
init_x = [C; D; D; INT; VS; VSf; GE];

%% Simulating the model
options = odeset('MaxStep', 0.1);
%[tinteg, full_states] = ode45(@(t,y) LaurensVestibularModelDeriv(t,y,u,params), [min(t), max(t)], init_x, options);
%[tinteg, full_states] = ode45(@(t,y) CanalDynamicsDeriv(t,y,u,params), [min(t), max(t)], C);%init_x);
%[tinteg, full_states] = ode45(@(t,y) LaurensVestibularModelDeriv(t,y,u,params), [min(t), max(t)], [C;D;INT;VS;VSf]);%init_x);

[tinteg, full_states] = ode45(@(t,y) SpineLabDualAdaptationDeriv(t,y,u,params), [min(t), max(t)], init_x, options);

% Get velocity estimate
nSteps = length( tinteg );
VE = zeros( 3, nSteps );
D = zeros( 3, nSteps );
AW = zeros( 2, nSteps );
for i=1:length(tinteg)
    
    % Get the fused canal signal
    Dinput = [full_states(i,4:6); full_states(i,7:9)]';
    [D(:,i), S(:,i), AW(:,i)] = CanalParticleWeighting( Dinput, [0,0,0], params.sigprior );
    
    u.D = D(:,i); %full_states(i,4:6)'; %
    u.dirVI = [0;0;0];
    u.VSf = full_states(i,16:18)';
    %u.VSf = full_states(i,13:15)';
    VE(:,i) = VelocityEstimate( tinteg(i), u, params );
end