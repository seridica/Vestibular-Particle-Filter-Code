% Laurens model
% details will follow later

clear all;
%close all;

dt = 0.1;   % Time step of the model
Fs = 1/dt;

%% Model parameters
Tcan = -1*[1 0 0;...
      -0.12 -0.7 0.7;...
      -0.12 -0.7 -0.7];% Vector normals to the semicircular canals (Canal orientation) row1 = lateral canal; row 2 = % anterior; row 3 = posterior
  
Tcan = [1 0 0;...
      0 1 0;...
      0 0 1];% Vector normals used to evaluate the response of canals oriented to the global system (mostly used for debugging)
  

grav = [0; 0; -9.81];  

% [B,A] = butter(3,(10/Fs/2),'low');
  
%% Motion profiles
mopo = 1;   % Evaluate motion responses recreating responses in Laurens and Angelaki paper
if mopo == 1
    % Parameters for Laurens and Angelaki EBR paper
    kv = 0.247;                      % Cupula gain (Laurens and Angelaki) - NOTE: this value was provided by Laurens by email in which he said 0.2 was incorrect
    Tc = 4;                          % Cupula time constant (Laurens and Angelaki)
    Tn=75;                           % Afferent adaptation time constant(George et al. 2011) 
%     Tr=15;                           % Adaptation recovery time constant (for now is just an arbitrary number)
    gd = 1;                          % Direct vestibular gain (Laurens and Angelaki)
    Tvs = 14;                        % Velocity storage time constant (Laurens and Angelaki)

    go = 2;                          % Direct visual gain (Laurens, Valko, etc., Laurens and Angelaki value is too high)
    ko = 0.6;                        % Indirect visual gain (Laurens and Angelaki)

    kf = 0.38;
    Ts = 1/0.65;
    d = 1;                           % Somatogravic feedback gain (not an actual parameter - used as on/off switch)

    sigC = [0.0]; %[0.1]*pi/180;             % Canal noise (Laurens and Angelaki)
    sigVS = [0.0]; %[0.2]*pi/180;            % Vestibular storage noise (Laurens and Angelaki seemed to high, halfed it)
    sigVI = [0.0]; %[0.1]*pi/180;            % Vision noise (Laurens and Angelaki seemed high, 1/10 of listed value)
    
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

    sigC = [0.1]*pi/180;       % Canal noise
    sigVS = [0.4]*pi/180;      % Vistublar storage noise
    sigVI = [1]*pi/180;        % Vision noise
    
    %%% Angular inputs
    % Pitch profile
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
elseif mopo == 3
    % Parameters for Laurens and Angelaki EBR paper
    kv = 0.247;                      % Cupula gain (Laurens and Angelaki) - NOTE: this value was provided by Laurens by email in which he said 0.2 was incorrect
    Tc = 4;                          % Cupula time constant (Laurens and Angelaki)
    Tn=75;                           % Afferent adaptation time constant(George et al. 2011) 
%     Tr=15;                           % Adaptation recovery time constant (for now is just an arbitrary number)
    gd = 1;                          % Direct vestibular gain (Laurens and Angelaki)
    Tvs = 14;                        % Velocity storage time constant (Laurens and Angelaki)

    go = 2;                          % Direct visual gain (Laurens, Valko, etc., Laurens and Angelaki value is too high)
    ko = 0.6;                        % Indirect visual gain (Laurens and Angelaki)

    kf = 0.38;
    Ts = 1/0.65;
    d = 1;                           % Somatogravic feedback gain (not an actual parameter - used as on/off switch)

    sigC = [0.0]; %[0.1]*pi/180;             % Canal noise (Laurens and Angelaki)
    sigVS = [0.0]; %[0.2]*pi/180;            % Vestibular storage noise (Laurens and Angelaki seemed to high, halfed it)
    sigVI = [0.0]; %[0.1]*pi/180;            % Vision noise (Laurens and Angelaki seemed high, 1/10 of listed value)
    
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
    
    % Time vector
    N = length(t);
else
end

%% Time vector
% Time vector
t = (0:length(angAcc)-1)*dt;
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

% Initializing output arrays
C_full = zeros(size(angAcc));
Cnn_full = zeros(size(angAcc));
C_full_A=zeros(size(angAcc));
Cnn_full_A=zeros(size(angAcc));   
VS_full = zeros(size(angAcc));
VE_full = zeros(size(angAcc));
GE_full = zeros(size(angAcc));
AE_full = zeros(size(angAcc));
GIA_full = zeros(size(angAcc));
Grav_full = zeros(size(angAcc));
Tori_full = zeros([3 size(angAcc)]);
angCAccT_full = zeros(size(angAcc));

INTnn_full = zeros(size(angAcc));
INT_full = zeros(size(angAcc));
VSli_full = zeros(size(angAcc));
VSvi_full = zeros(size(angAcc));
VEi_full = zeros(size(angAcc));
VEli_full = zeros(size(angAcc));
VEvi_full = zeros(size(angAcc));

%% Simulating the model - difference approach

angCAcc = Tcan*angAcc;           % Transform angular rotation to canal coordinates

for rr = 1                       % Loop for running several iterations if desired
    C = zeros(3,1);              % Canal output signal   
    DA = zeros(3,1);              % neural adaptation variable
    DAnn = zeros(3,1);
    Cnn = zeros(3,1);            % Ideal canal output signal - no noise
    VS = zeros(3,1);             % Velocity storage
    VE = zeros(3,1);             % Velocity estimate - angular

    INTnn = zeros(3,1);          % Ideal integral - no noise from anything
    INT = zeros(3,1);            % Integral output - with noise (canals and VS)
    VSli = zeros(3,1);           % Velocity storage - no vision
    VEi = zeros(3,1);            % Velocity estimate - noisy integrator
    VEli = zeros(3,1);           % Velocity estimate - leaky integrator
    
    VSnoise_prev = zeros(3,1);   % Noise from previous step 
    
    GE = rotx(angPos(1,1)) * roty(angPos(2,1)) * rotz(angPos(3,1)) * grav;
    GE_full(:,1,rr) = GE;
    Grav_full(:,1,rr) = GE;
    Tori_full(:,:,1,rr) = rotx(angPos(1,1)) * roty(angPos(2,1)) * rotz(angPos(3,1));
    for k = 2:length(t)
        % Noise values
        Cnoise = sigC*randn(3,1);
        VSnoise = sigVS*randn(3,1);
        dVSnoise = VSnoise - VSnoise_prev;
        VSnoise_prev = VSnoise;
        VInoise = sigVI*randn(3,1);
        
        % Gravity input signals
        Tori = rotx(angPos(1,k)) * roty(angPos(2,k)) * rotz(angPos(3,k));
        Tori_full(:,:,k) = Tori;
        gravT = Tori * grav;
        Grav_full(:,k) = gravT;
        
        tranAccT = Tori * tranAcc(:,k);
        GIA = gravT - tranAccT;
        GIA_full(:,k,rr) = GIA;
        
        % Cupula equations
        angAccT = Tori*angAcc(:,k);
        angCAccT = Tcan*angAccT;
        angCAccT_full(:,k) = angCAccT;
        Cnn = Cnn*(1-dt/Tc) + angCAccT*dt;
        C = C*(1-dt/Tc) + angCAccT*dt + Cnoise;
        %Cnn = Cnn*exp(-dt/Tc) + angCAccT*dt;
        %C = C*(exp(-dt/Tc)) + angCAccT*dt + Cnoise;
        % implementing adaptation 
       
        
% % % %         if C>0.02           
% % % %         DA=DA+C/Tn*dt;
% % % %         elseif 0<C<0.02
% % % %             DA=DA-DA/Tr*dt;
% % % %         end;
    
        
% % % % % % % % %         D=D*(1-exp(-dt/Tn));
% % % % %         C=D*C;
        C_full(:,k) = C;   % Cupula output signal
%         DA=C_full(:,1);
%         for ii=2:k
%              DA=DA*exp(-dt/Tn)+C_full(:,ii)-C_full(:,ii-1);
%         end
        DA = DA*(1-dt/Tn) + (C_full(:,k)-C_full(:,k-1));
        C_full_A(:,k)=DA; % Afferent output (with adaptation) 
        
        % feeding other parts of the model with the presence of afferent
        % adaptation
        Cnn_full(:,k) = Cnn;
        DAnn = DAnn*(1-dt/Tn) + (Cnn_full(:,k)-Cnn_full(:,k-1));
%          DAnn=Cnn_full(:,1);
%         for ii=2:k
%              DAnn=DAnn*exp(-dt/Tn)+Cnn_full(:,ii)-Cnn_full(:,ii-1);
%         end
        Cnn_full_A(:,k)=DAnn;
        % Integrator equations (ideal endolymph and noisy endolymph positions)
        INTnn = INTnn + kv*DAnn*dt; 
        INTnn_full(:,k,rr) = INTnn;
        INT = INT + (kv*DA)*dt + dVSnoise;        % noise is on the output of the integral
        INT_full(:,k,rr) = INT;        

        % Retinal slip
%         rSL = INTnn - VSli;
% Retinal slip according to patricks email:
         rSL = (angVel(:,k)+VInoise-VSli-DA)/(1+go);  
        % Indirect visual pathway
        indVI = (rSL + VInoise)*ko;
        % Direct visual pathway
        dirVI = (rSL + VInoise)*go;

        % Velocity storage
        %VS = VS*exp(-dt/Tvs) + (kv*DA + indVI + kf*Tcan*cross(GIA/norm(GIA),GE/norm(GE)))*dt + dVSnoise;      % noise is on the output of the integral, therefore it is dVSnoise and not VSnoise
        VS = VS*(1-dt/Tvs) + (kv*DA + indVI + kf*Tcan*cross(GIA/norm(GIA),GE/norm(GE)))*dt + dVSnoise;      % noise is on the output of the integral, therefore it is dVSnoise and not VSnoise
        VS_full(:,k,rr) = VS; % Velocity storage output signal

        % Total rotational estimate
        VEprev = VE;
        VE = VS + gd*DA + dirVI;
        VE_full(:,k,rr) = VE;     % Complete velocity estimate signal
        % Gravity and acceleration estimates
        GE = (cross(GE,inv(Tcan)*VE))*dt + GE*(1-dt/Ts) + 1/Ts*GIA*dt; %GE*exp(-dt/Ts) + 1/Ts*GIA*dt;
        %GE = (cross(GE,inv(Tcan)*VE))*dt + GE*exp(-dt/Ts) + 1/Ts*GIA*dt;        
        GE_full(:,k,rr) = GE;
        AE = GIA - GE;       
        AE_full(:,k,rr) = AE;
        
        %%%%%%%%%%% Extra stuff for plotting Laurens and Angelaki paper 
        % Retinal slip - no gravitational contribution
        rSLvi = INTnn - VSvi;
        % Indirect visual pathway
        indVIvi = (rSLvi + VInoise)*ko ;
        % Direct visual pathway
        dirVIvi = (rSLvi + VInoise)*go;

        % Velocity estimate - integral only
        VEi = INT + gd*DA;                                       
        VEi_full(:,k,rr) = VEi;

        % Velocity storage and estimate - leaky integrator only - without vision and without gravity
        %VSli = VSli*exp(-dt/Tvs) + (kv*DA)*dt + dVSnoise;        % noise is on the output of the integral
        VSli = VSli*(1-dt/Tvs) + (kv*DA)*dt + dVSnoise;        % noise is on the output of the integral
        VSli_full(:,k,rr) = VSli;
        VEli = VSli + gd*DA;
        VEli_full(:,k,rr) = VEli;
        
        % Velocity storage and estimate - leaky integrator and vision - without gravity
        %VSvi = VSvi*exp(-dt/Tvs) + (kv*DA + indVIvi)*dt + dVSnoise;% noise is on the output of the integral
        VSvi = VSvi*(1-dt/Tvs) + (kv*DA + indVIvi)*dt + dVSnoise;% noise is on the output of the integral
        VSvi_full(:,k,rr) = VSvi;
        VEvi = VSvi + gd*DA+ dirVIvi;
        VEvi_full(:,k,rr) = VEvi;
    end
end

%% Plotting figure 1 from Laurens and Angelaki
if mopo == 1 || mopo == 3
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
    
    hC_full = inv(Tcan)*C_full*180/pi;
    hC_full_A = inv(Tcan)*C_full_A*180/pi;
    hVEi_full = inv(Tcan)*VEi_full*180/pi;
    hVEli_full = inv(Tcan)*VEli_full*180/pi;
    hVEvi_full = inv(Tcan)*VEvi_full*180/pi;
    hVE_full = inv(Tcan)*VE_full*180/pi;
    
%     INT_full = trnfrmTimsig(Tori_full,INT_full,1);
%     INTnn_full = trnfrmTimsig(Tori_full,INTnn_full,1);
%     VSli_full = trnfrmTimsig(Tori_full,VSli_full,1);
%     VSvi_full = trnfrmTimsig(Tori_full,VSvi_full,1);
%     VS_full = trnfrmTimsig(Tori_full,VS_full,1);
    
    fig = figure;
    angVel = angVel*180/pi;
    liu = 40; lil = -6;
    liu2 = 0.6; lil2 = -0.4;
    xli = Tmo/2;
 
        
    %%% Head in space
    % Canals
    set(fig,'Color','w');
    subplot(2,5,1);
    plot(t,angVel,'b','LineStyle',':'); hold on;
    plot(t,hC_full,'k');
    plot(t,hC_full_A, 'r');
    box off;
%     xlim([0 xli]);
%     ylim([lil liu]);
    title('Canals');
    ylabel('Head in space [\circ/s]');
    % Noisy integrator
    subplot(2,5,2);
    plot(t,angVel,'b','LineStyle',':'); hold on;
    plot(t,hVEi_full,'k');
    box off;
%     xlim([0 xli]);
%     ylim([lil liu]);
    title('Noisy integrator');
    % Leaky integrator
    subplot(2,5,3);
    plot(t,angVel,'b','LineStyle',':'); hold on;
    plot(t,hVEli_full,'k');
    box off;
%     xlim([0 xli]);
%     ylim([lil liu]);
    title('Leaky integrator'); 
    % Vision
    subplot(2,5,4);
    plot(t,angVel,'b','LineStyle',':'); hold on;
    plot(t,hVEvi_full,'k');
    box off;
%     xlim([0 xli]);
%     ylim([lil liu]);
    title('Vision');
    % Somatogravic feedback
    subplot(2,5,5);
    plot(t,angVel,'b','LineStyle',':'); hold on;
    plot(t,hVE_full,'k');
    box off;
%     xlim([0 xli]);
%     ylim([lil liu]);
    title('Somatogravic');

    %%% Endolymph in space
    % Canals
    subplot(2,5,6);
    plot(t,INTnn_full,'b','LineStyle',':');
    box off;
%     xlim([0 xli]);
%     ylim([lil2 liu2]);
    xlabel('Time [sec]');
    % Noisy integrator
    subplot(2,5,7);
    plot(t,INTnn_full,'b','LineStyle',':'); hold on;
    plot(t,INT_full,'k');
    box off;
%     xlim([0 xli]);
%     ylim([lil2 liu2]);
    xlabel('Time [sec]');
    % Leaky integrator
    subplot(2,5,8);
    plot(t,INTnn_full,'b','LineStyle',':'); hold on;
    plot(t,VSli_full,'k');
    box off;
%     xlim([0 xli]);
%     ylim([lil2 liu2]);
    xlabel('Time [sec]');
    % Vision
    subplot(2,5,9);
    plot(t,INTnn_full,'b','LineStyle',':'); hold on;
    plot(t,VSvi_full,'k');
    box off;
%     xlim([0 xli]);
%     ylim([lil2 liu2]);
    xlabel('Time [sec]');
    % Somatogravic
    subplot(2,5,10);
    plot(t,INTnn_full,'b','LineStyle',':'); hold on;
    plot(t,VS_full,'k');
    box off;
%     xlim([0 xli]);
%     ylim([lil2 liu2]);
    xlabel('Time [sec]');
end

% Plotting here is just for viewing, doesn't match any papers
% Gravitational estimate
figure;
plot(t,GE_full,'Linewidth',2); hold on;
ylim([-12 12]);
plot(t,angVel,'k');
plot(t,tranAcc,'r');
ylabel('GE');
legend({'x','y','z'});

% Estimated angle
ang = atan(GE_full(1,:)./GE_full(3,:))*180/pi;
figure;
plot(t,ang); hold on;
plot(t,angVel,'k');
ylabel('Estimated angle');

% Acceleration esimates
figure;
plot(t,AE_full,'Linewidth',2); hold on;
ylabel('AE');
legend({'x','y','z'});
