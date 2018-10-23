%%%% %%%%%%%%%%%%% Laurens vestibular model
%
% Name: vestibular internal model based on Laurens et al. 2013, Nature
% Neuroscience paper
% Authors: PAF, CJD, NK, JSB
% Purpose: simulate model desribed in Nature Neuroscience paper
% Date: 13-01-2017
%
%%%%  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all; 
close all;

%spacelib

mopo=2;
dt = 0.1;   % Time step of the model
Fs = 1/dt;

%% Model parameters
% Vector normals to the semicircular canals (Canal orientation) 
% row 1 = % anterior; row 2 = posterior; row 3 = lateral canal;
Tcan = -1*[-0.12 -0.7 0.7;...
          -0.12 -0.7 -0.7;...
          1      0    0];

grav = [0; 0; -9.81];  

% [B,A] = butter(3,(10/Fs/2),'low');
  
%% Motion profiles - angular inputs (no translation for Laurens model)

% Pitch rotation of the Head in Chair Space
stpT   = 1.4;
hldT   = 32-stpT;
stp    = -1+2/(stpT/dt):2/(stpT/dt):1-2/(stpT/dt);
stp2   = [0    0.1429    0.2857    0.4286    0.5714    0.7143    0.8571];
pitPos = 10*[0 0 stp2 ones(1,10/dt) fliplr(stp) -ones(1,round(hldT/dt)+1) stp ones(1,round(hldT/dt)+1) fliplr(stp) -ones(1,round(hldT/dt)+1) stp ones(1,round(hldT/dt)+1)]; 
pitVel = [0, diff(pitPos)]/dt;
    
% Yaw rotation of the Chair in Lab Space
lengPit = length(pitVel);
yawVel  = [zeros(1,8/dt) linspace(0,45*pi/180,20) 45*pi/180*ones(1,lengPit-length(zeros(1,8/dt))-length(linspace(0,45*pi/180,20)))]; 
tt      = (0:length(yawVel)-1)*dt;
yawPos  = cumtrapz(tt,yawVel);
    
% 3D vectors for Yaw (in Lab Space) and Pitch (in Chair Space) positions in Rad
angPos_yaw = [zeros(1,lengPit); zeros(1,lengPit); yawPos]; % position matrice of yaw movement in universal coordinate
angPos_pit = [zeros(1,lengPit); pitPos; zeros(1,lengPit)]*pi/180; % position matrice of pitch movement in head coordinate

% 3D vectors for Yaw (in Lab Space) and Pitch (in Chair Space) velocities
angVel_yaw = [zeros(1,lengPit); zeros(1,lengPit); yawVel];
angVel_pit = [zeros(1,lengPit); pitVel*pi/180; zeros(1,lengPit)];

%% Rotation matrices and 3D tranformations - full input arrays

% Yaw velocity in Head space - Vel instead of Pos to avoid discontinuities
for kk = 2:length(tt)
    % Generating the rotation matrice from chair to head coordinates for the pitch movement (angPos_pit)    
    rotm_L2H(:,:,kk) = cardator(angPos_pit(:,kk)',X, Y, Z); 
    % Yaw velocity 
    YawinChair(:,kk) = angVel_yaw(:,kk)' * rotm_L2H(:,:,kk);
    %Orientation of gravity in the head
    gravT(:,kk) = (rotm_L2H(3,1:3,kk))' * -9.81; %Determine the direction of + gravitational field
end
gravT(:,1) = [0 0 -9.81];

% Adding Pitch velocity of the Head to Yaw motion in Chair space
angVel_Chair = YawinChair + angVel_pit;

% Head and canal velocity in Lab space - Vel instead of Pos to avoid discontinuities
for kk=2:length(tt)
   rotm_L2C = cardator(angPos_yaw(:,kk)',X,Y,Z); % generating the rotation matrice from universal to chair coordinate for the yaw movement
   angVel_Lab(:,kk) = angVel_Chair(:,kk)' * inv(rotm_L2C);
   angVel_Can(:,kk) = angVel_Chair(:,kk)' * Tcan;
end

% Generating Angular position and angular velocity variables - in Space Chair 
t = (0:length(angVel_Lab)-1)*dt;
angPos_Lab   = cumtrapz(t,angVel_Lab,2);
angPos_Chair = cumtrapz(t,angVel_Chair,2);
angPos_Can   = cumtrapz(t,angVel_Can,2);
angAcc_Lab   = [zeros(3,1), diff(angVel_Lab,[],2)]/dt;
angAcc_Chair = [zeros(3,1), diff(angVel_Chair,[],2)]/dt;
angAcc_Can   = [zeros(3,1), diff(angVel_Can,[],2)]/dt;

% Translations
tranAcc = zeros(size(angAcc_Chair));

%% Model parameters - all parameters mimic Nature Neuroscience paper
kv = 0.2;                  % Cupula gain
Tc = 4;                    % Cupula time constant 
gd = 1;                    % Direct vestibular gain 
Tvs = 10;                  % Velocity storage time constant 

go = 0;                    % Direct visual gain = set to zero for this sim
ko = 0.0;                  % Indirect visual gain = set to zero for this sim

% kf = 0.05;                 % Velocity feedback gain (in nature supplementary materials)
kf = 0.5;                 % Velocity feedback gain (seems this one gives a closer resemblance to nature figures)
Ts = 0.5;                  % Somatogravic feedback time constant
d = 1;                     % Somatogravic feedback gain

sigC = [0]*pi/180;         % Canal noise
sigVS = [0]*pi/180;        % Vistublar storage noise
sigVI = [0]*pi/180;        % Vision noise
    
%% Time vector
% Time vector
t = (0:length(angAcc_Chair)-1)*dt;
N = length(t);

%% Initializing the variables
% Initializing initial conditions
C = zeros(3,1);              % Canal output signal   
Cnn = zeros(3,1);            % Ideal canal output signal - no noise
VS = zeros(3,1);             % Velocity storage
VE = zeros(3,1);             % Velocity estimate - angular
GE = zeros(3,1);             % Gravity estimate
AE = zeros(3,1);             % Acceleration estimate
GIA = zeros(3,1);            % Gravitoinertial acceleration

INTnn = zeros(3,1);          % Ideal integral - no noise from anything
INT = zeros(3,1);            % Integral output - with noise (canals and VS)
VSli = zeros(3,1);           % Velocity storage - leaky integrator
VSvi = zeros(3,1);           % Velocity storage - vision
VEi = zeros(3,1);            % Velocity estimate - noisy integrator
VEli = zeros(3,1);           % Velocity estimate - leaky integrator
VEvi = zeros(3,1);           % Velocity estimate - vision
GEli = zeros(3,1);           % Gravity estimate - leaky integrator only (no somatogravic or rotational feedback)
AEli = zeros(3,1);           % Acceleration estimate - leaky integrator only (no somatogravic or rotational feedback)
GEsf = zeros(3,1);           % Gravity estimate - leaky integrator with somatogravic feedback (no rotational feedback)
AEsf = zeros(3,1);           % Acceleration estimate - leaky integrator with somatogravic feedback  (no rotational feedback)

% Initializing output arrays - see parameter definitions above
C_full = zeros(size(angAcc_Chair));
Cnn_full = zeros(size(angAcc_Chair));
C_full_A=zeros(size(angAcc_Chair));
Cnn_full_A=zeros(size(angAcc_Chair));   
VS_full = zeros(size(angAcc_Chair));
VE_full = zeros(size(angAcc_Chair));
GE_full = zeros(size(angAcc_Chair));
AE_full = zeros(size(angAcc_Chair));
GIA_full = zeros(size(angAcc_Chair));

INTnn_full = zeros(size(angAcc_Chair));
INT_full = zeros(size(angAcc_Chair));
VSli_full = zeros(size(angAcc_Chair));
VSvi_full = zeros(size(angAcc_Chair));
VEi_full = zeros(size(angAcc_Chair));
VEli_full = zeros(size(angAcc_Chair));
VEvi_full = zeros(size(angAcc_Chair));
GEli_full = zeros(size(angAcc_Chair));
AEli_full = zeros(size(angAcc_Chair));
GEsf_full = zeros(size(angAcc_Chair));
AEsf_full = zeros(size(angAcc_Chair));

%% Simulating the model - difference approach

for rr = 1                       % Loop for running several iterations if desired
    C = zeros(3,1);              % Canal output signal   
    Cnn = zeros(3,1);            % Ideal canal output signal - no noise
    VS = zeros(3,1);             % Velocity storage
    VE = zeros(3,1);             % Velocity estimate - angular

    INTnn = zeros(3,1);          % Ideal integral - no noise from anything
    INT = zeros(3,1);            % Integral output - with noise (canals and VS)
    VSli = zeros(3,1);           % Velocity storage - no vision
    VEi = zeros(3,1);            % Velocity estimate - noisy integrator
    VEli = zeros(3,1);           % Velocity estimate - leaky integrator
    
    VSnoise_prev = zeros(3,1);   % Noise from previous step 
    
    GE = gravT(:,1);
    GEli = GE;
    GEsf = GE;
    GE_full(:,1,rr) = GE;
    GEli_full(:,1,rr) = GE; 
    GEsf_full(:,1,rr) = GE;
    for k = 2:length(t)
        % Noise values
        Cnoise = sigC*randn(3,1);
        VSnoise = sigVS*randn(3,1);
        dVSnoise = VSnoise - VSnoise_prev;
        VSnoise_prev = VSnoise;
        VInoise = sigVI*randn(3,1);
        
        % Gravity and translation input signals - calculate GIA
        GIA = gravT(:,k) - tranAcc(:,k);
        GIA_full(:,k,rr) = GIA;
        
        % Cupula equations
        angAcc = angAcc_Chair(:,k);
        Cnn = Cnn*exp(-dt/Tc) + angAcc*dt; 
        C = C*(exp(-dt/Tc)) + angAcc*dt + Cnoise;
        C_full(:,k,rr) = C;   % Cupula output signal

        % Integrator equations (ideal endolymph and noisy endolymph positions)
        INTnn = INTnn + kv*Cnn*dt; 
        INTnn_full(:,k,rr) = INTnn;
        INT = INT + (kv*C)*dt + dVSnoise;        % noise is on the output of the integral
        INT_full(:,k,rr) = INT;        

        % Retinal slip - equations obtained from Jean Laurens - email 30-10-2016 sent to PAF
        rSL = (angVel_Chair(:,k) + VInoise - VS - C)/(1+go);  
        % Indirect visual pathway
        indVI = (rSL)*ko;
        % Direct visual pathway
        dirVI = (rSL)*go;

        % Velocity storage
        VS = VS*exp(-dt/Tvs) + (kv*C + indVI + kf*cross(GIA/norm(GIA),GE/norm(GE)))*dt + dVSnoise;      % noise is on the output of the integral, therefore it is dVSnoise and not VSnoise
        VS_full(:,k,rr) = VS; % Velocity storage output signal

        % Total rotational estimate
        VEprev = VE;
        VE = VS + gd*C + dirVI;
        VE_full(:,k,rr) = VE;     % Complete velocity estimate signal
        
        % Gravity and acceleration estimates    %%%%%%% MAY BE DIFFICULT FOR SIMULINK
        Amat = [-1/Ts VE(3) -VE(2); -VE(3) -1/Ts VE(1); VE(2) -VE(1) -1/Ts];
        GE = expm(Amat*dt)*GE + inv(Amat)*(expm(Amat*dt)-eye(3))*1/Ts*GIA;
        GE_full(:,k,rr) = GE;
        AE = GE-GIA;       
        AE_full(:,k,rr) = AE;
        
        %%%%%%%%%%%%% WARNING WARNING WARNING WARNING %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Extra stuff for plotting Laurens and Angelaki EBR paper 
        % UNLESS SPECIFICALLY PLOTTING EBR GRAPHS, THEN DON'T USE THESE VARIABLES!!!!!!!!! 
        % Retinal slip - no gravitational contribution
        rSLvi = (angVel_Chair(:,k)+VInoise-VSli-C)/(1+go);  
        % Indirect visual pathway
        indVIvi = (rSLvi)*ko;
        % Direct visual pathway
        dirVIvi = (rSLvi)*go;

        % Velocity estimate - integral only    
        VEi = INT + gd*C;                                       
        VEi_full(:,k,rr) = VEi;

        % Velocity storage and estimate - leaky integrator only - without vision and without gravity
        VSli = VSli*exp(-dt/Tvs) + (kv*C)*dt + dVSnoise;        % noise is on the output of the integral
        VSli_full(:,k,rr) = VSli;
        VEli = VSli + gd*C;
        VEli_full(:,k,rr) = VEli;
        
        % Velocity storage and estimate - leaky integrator and vision - without gravity
        VSvi = VSvi*exp(-dt/Tvs) + (kv*C + indVIvi)*dt + dVSnoise;% noise is on the output of the integral
        VSvi_full(:,k,rr) = VSvi;
        VEvi = VSvi + gd*C + dirVIvi;
        VEvi_full(:,k,rr) = VEvi;
        
        % Gravity estimate without vision, somatogravic feedback and rotational feedback - second part of Figure 2 EBR paper
        Amat = [0 VEli(3) -VEli(2); -VEli(3) 0 VEli(1); VEli(2) -VEli(1) 0];
        GEli = expm(Amat*dt)*GEli;
        GEli_full(:,k,rr) = GEli;
        AEli = GIA - GEli;       
        AEli_full(:,k,rr) = AEli;
        
        % Gravity estimate without vision and rotational feedback - third part of Figure 2 EBR paper
        Amat = [-1/Ts VEli(3) -VEli(2); -VEli(3) -1/Ts VEli(1); VEli(2) -VEli(1) -1/Ts];
        GEsf = expm(Amat*dt)*GEsf + inv(Amat)*(expm(Amat*dt)-eye(3))*1/Ts*GIA;
        GEsf_full(:,k,rr) = GEsf;
        AEsf = GIA - GEsf;       
        AEsf_full(:,k,rr) = AEsf;
    end
end

%% Plotting figure 1 from Laurens and Angelaki
C_full = mean(C_full,3);
VEi_full = mean(VEi_full,3);
VEli_full = mean(VEli_full,3);
VEvi_full = mean(VEvi_full,3);
VE_full = mean(VE_full,3);

INT_full = mean(INT_full,3);
INTnn_full = mean(INTnn_full,3);
VSli_full = mean(VSli_full,3);
VSvi_full = mean(VSvi_full,3);
VS_full = mean(VS_full,3);

GIA_full = mean(GIA_full,3);
GE_full = mean(GE_full,3);
GEli_full = mean(GEli_full,3);
GEsf_full = mean(GEsf_full,3);
AE_full = mean(AE_full,3);
AEli_full = mean(AEli_full,3);
AEsf_full = mean(AEsf_full,3);

% Extra signals for EBR paper, these need updating to be used for making
% the plot below.  Best not to use them for now.

% hC_full = inv(Tcan)*C_full*180/pi;
% hVEi_full = inv(Tcan)*VEi_full*180/pi;
% hVEli_full = inv(Tcan)*VEli_full*180/pi;
% hVEvi_full = inv(Tcan)*VEvi_full*180/pi;
% hVE_full = inv(Tcan)*VE_full*180/pi;

% INT_full = trnfrmTimsig(Tori_full,INT_full,1);
% INTnn_full = trnfrmTimsig(Tori_full,INTnn_full,1);
% VSli_full = trnfrmTimsig(Tori_full,VSli_full,1);
% VSvi_full = trnfrmTimsig(Tori_full,VSvi_full,1);
% VS_full = trnfrmTimsig(Tori_full,VS_full,1);
    
%% Plotting the Nature Neuroscience Figure 2
fig = figure;
set(fig,'Color','w','Units','Normalized','Position',[0.1 0.1 0.6 0.5]);
% Plot pitch position of chair
subplot(5,1,1);
plot([0 150],[0 0],'k--'); hold on;
plot(t,angPos_Chair(2,:)*180/pi,'g');
ylim([-10 10]);
box off;
ylabel('Pitch (\circ)');
% Plot pitch velocity of the chair
subplot(5,1,2);
plot([0 150],[0 0],'k--'); hold on;
plot(t,angVel_Chair(2,:)*180/pi,'g');
ylim([-20 20]);
box off;
ylabel('Pitch vel (\circ/s)');
% Plot the yaw velocity of the chair, canals and rotation estimate
subplot(5,1,3);
plot([0 150],[0 0],'k--'); hold on;    
plot(t,angVel_Chair(3,:)*180/pi,'b','LineStyle',':');
plot(t,C_full(3,:)*180/pi,'k');
plot(t,VE_full(3,:)*180/pi,'b');
ylim([0 45]);
box off;
ylabel('Yaw vel (\circ/s)');
% Plot the roll velocity estimate
subplot(5,1,4);
plot([0 150],[0 0],'k--'); hold on;
plot(t,angVel_Chair(1,:)*180/pi,'c:');
plot(t,VE_full(1,:)*180/pi,'c');
ylim([-15 15]);
box off;
ylabel('Roll vel (\circ/s)');
subplot(5,1,5);
plot([0 150],[0 0],'k--'); hold on;  
plot(t,tranAcc(2,:),'r:');
plot(t,AE_full(2,:),'r');
% plot(t,AE_full([1 3],:));
ylim([-1.2 1.2]);
box off;
ylabel('Accel est [m/s^2]');
   
%% These plotting scripts are for replicating the EBR papers - ignore if not needed
mopo = 0;
if mopo == 1
    fig = figure;
    angVel = angVel*180/pi;
    liu = 40; lil = -20;
    liu2 = 0.6; lil2 = -0.4;
    xli = Tmo;
        
    %%% Head in space
    % Canals
    set(fig,'Color','w','Units','Normalized','Position',[0.1 0.1 0.6 0.5]);
    subplot(2,5,1);
    plot(t,angVel,'b','LineStyle',':'); hold on;
    plot(t,hC_full,'k');
    box off;
    xlim([0 xli]);
    ylim([lil liu]);
    title('Canals');
    ylabel('Head in space [\circ/s]');
    % Noisy integrator
    subplot(2,5,2);
    plot(t,angVel,'b','LineStyle',':'); hold on;
    plot(t,hVEi_full,'k');
    box off;
    xlim([0 xli]);
    ylim([lil liu]);
    title('Noisy integrator');
    % Leaky integrator
    subplot(2,5,3);
    plot(t,angVel,'b','LineStyle',':'); hold on;
    plot(t,hVEli_full,'k');
    box off;
    xlim([0 xli]);
    ylim([lil liu]);
    title('Leaky integrator'); 
    % Vision
    subplot(2,5,4);
    plot(t,angVel,'b','LineStyle',':'); hold on;
    plot(t,hVEvi_full,'k');
    box off;
    xlim([0 xli]);
    ylim([lil liu]);
    title('Vision');
    % Somatogravic feedback
    subplot(2,5,5);
    plot(t,angVel,'b','LineStyle',':'); hold on;
    plot(t,hVE_full,'k');
    box off;
    xlim([0 xli]);
    ylim([lil liu]);
    title({'Somatogravic and';'rotation feedback'});

    %%% Endolymph in space
    % Canals
    subplot(2,5,6);
    plot(t,INTnn_full,'b','LineStyle',':');
    box off;
    xlim([0 xli]);
    ylim([lil2 liu2]);
    xlabel('Time [sec]');
    ylabel('Endolymph [\circ/s]');
    % Noisy integrator
    subplot(2,5,7);
    plot(t,INTnn_full,'b','LineStyle',':'); hold on;
    plot(t,INT_full,'k');
    box off;
    xlim([0 xli]);
    ylim([lil2 liu2]);
    xlabel('Time [sec]');
    % Leaky integrator
    subplot(2,5,8);
    plot(t,INTnn_full,'b','LineStyle',':'); hold on;
    plot(t,VSli_full,'k');
    box off;
    xlim([0 xli]);
    ylim([lil2 liu2]);
    xlabel('Time [sec]');
    % Vision
    subplot(2,5,9);
    plot(t,INTnn_full,'b','LineStyle',':'); hold on;
    plot(t,VSvi_full,'k');
    box off;
    xlim([0 xli]);
    ylim([lil2 liu2]);
    xlabel('Time [sec]');
    % Somatogravic
    subplot(2,5,10);
    plot(t,INTnn_full,'b','LineStyle',':'); hold on;
    plot(t,VS_full,'k');
    box off;
    xlim([0 xli]);
    ylim([lil2 liu2]);
    xlabel('Time [sec]');
    
elseif mopo == 2
    fig = figure;
    set(fig,'Color','w','Units','Normalized','Position',[0.1 0.1 0.6 0.5]);
    angPos = angPos*180/pi;
    liu = 3; lil = -2;
    liu2 = 15; lil2 = -5;
    xli = Tmo;
    
    % Head orientation 
    subplot(3,4,1);
    plot(t,-angPos);
    ylabel('Tilt r.e. gravity [\circ]');
    box off;
    xlim([0 xli]);
    ylim([lil2 liu2]);
    title('Input');
    % Acceleration of head
    subplot(3,4,5);
    plot(t,-tranAcc);    
    box off;
    xlim([0 xli]);
    ylim([lil liu]);
    ylabel({'Acceleration';'(sign inverted)'});
    % GIA
    subplot(3,4,9);
    plot(t,GIA_full);
    box off;
    xlim([0 xli]);
    ylim([lil liu]);
    xlabel('Time [sec]');
    ylabel('GIA [m/s^2]');

    % Noisy tilt estimate
    subplot(3,4,2);
    ang = -atan(GEli_full(2,:)./GEli_full(3,:))*180/pi;    
    plot(t,ang);    
    box off;
    xlim([0 xli]);
    ylim([lil2 liu2]);
    title('Noisy estimate');
    % Acceleration of head
    subplot(3,4,6);
    plot(t,AEli_full);    
    box off;
    xlim([0 xli]);
    ylim([lil liu]);
    % GIA
    subplot(3,4,10);
    plot(t,GIA_full);
    box off;
    xlim([0 xli]);
    ylim([lil liu]);
    xlabel('Time [sec]');
    
    % Somatogravic feedback
    subplot(3,4,3);
    ang = -atan(GEsf_full(2,:)./GEsf_full(3,:))*180/pi;    
    plot(t,ang);    
    box off;
    xlim([0 xli]);
    ylim([lil2 liu2]);
    title({'Somatogravic'; 'feedback'});
    % Acceleration of head
    subplot(3,4,7);
    plot(t,AEsf_full);    
    box off;
    xlim([0 xli]);
    ylim([lil liu]);
    % GIA
    subplot(3,4,11);
    plot(t,GIA_full);
    box off;
    xlim([0 xli]);
    ylim([lil liu]);
    xlabel('Time [sec]');
    
    % Somatogravic + rotation feedback
    subplot(3,4,4);
    ang = -atan(GE_full(2,:)./GE_full(3,:))*180/pi;    
    plot(t,ang);    
    box off;
    xlim([0 xli]);
    ylim([lil2 liu2]);
    title({'Somatogravic and';'rotation feedback'});
    % Acceleration of head
    subplot(3,4,8);
    plot(t,AE_full);    
    box off;
    xlim([0 xli]);
    ylim([lil liu]);
    legend({'x','y','z'},'Location','NorthEast');
    legend boxoff;
    % GIA
    subplot(3,4,12);
    plot(t,GIA_full);
    box off;
    xlim([0 xli]);
    ylim([lil liu]);
    xlabel('Time [sec]');    
elseif mopo == 3
    fig = figure;
    set(fig,'Color','w','Units','Normalized','Position',[0.1 0.1 0.6 0.5]);
    angPos = angPos*180/pi;
    angVel = angVel*180/pi;
    liu = 5; lil = -5;
    liu2 = 95; lil2 = 80;
    xli = Tmo;
    
    % Input velocity
    subplot(2,4,1);
    plot(t,angVel);
    box off;
    xlim([0 xli]);
    ylim([lil liu]);
    ylabel('Velocity [\circ/s]');   
    title('Input');
    % No feedback velocity estimate
    subplot(2,4,2);
    plot(t,hVEli_full);
    box off;
    xlim([0 xli]);
    ylim([lil liu]);
    title('No feedback');
    % Somatogravic feedback
    subplot(2,4,3);
    plot(t,hVEli_full);
    box off;
    xlim([0 xli]);
    ylim([lil liu]);
    title({'Somatogravic'; 'feedback'});
    % Somatogravic and rotation feedback
    subplot(2,4,4);
    plot(t,hVE_full);
    box off;
    xlim([0 xli]);
    ylim([lil liu]);
    title({'Somatogravic and';'rotation feedback'});
    legend({'x','y','z'},'Location','NorthEast');
    legend boxoff;
    
    % Input orientation
    subplot(2,4,5);
    plot(t,angPos);
    box off;
    xlim([0 xli]);
    ylim([lil2 liu2]);
    ylabel('Tilt [\circ]');    
    xlabel('Time [sec]');
    % No feedback tilt estimate
    subplot(2,4,6);
    ang = -atan(GEli_full(1,:)./GEli_full(3,:))*180/pi;
    ind = find(ang < 0);
    ang(ind) = ang(ind)+180;    
    plot(t,(ang));    
    box off;
    xlim([0 xli]);
    ylim([lil2 liu2]);
    xlabel('Time [sec]');
    % Somatogravic feedback
    subplot(2,4,7);
    ang = -atan(GEsf_full(1,:)./GEsf_full(3,:))*180/pi;
    ind = find(ang < 0);
    ang(ind) = ang(ind)+180;    
    plot(t,(ang));
    box off;
    xlim([0 xli]);
    ylim([lil2 liu2]);    
    xlabel('Time [sec]');    
    % Somatogravic + rotation feedback
    subplot(2,4,8);
    ang = -atan(GE_full(1,:)./GE_full(3,:))*180/pi;
    ind = find(ang < 0);
    ang(ind) = ang(ind)+180;    
    plot(t,(ang));
    box off;
    xlim([0 xli]);
    ylim([lil2 liu2]);  
    xlabel('Time [sec]');    
end
   