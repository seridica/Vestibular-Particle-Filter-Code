function [t, angAcc] = MotionProfile( mopo, dt, dur )
    Fs = 1/dt;
    switch mopo
        %% Constant Angular Acceleration-Deceleration
        case 1
            angVMag = [0; 900; 0];
            Tmo = 220;
            t =-50:dt:Tmo;
            angPro = zeros(1,length(t));
            angPro(t>1&t<=76) = linspace(0,1,length(angPro(t>1&t<=76)));
            angPro(t>76&t<116)= 1;
            angPro(t>=116&t<191) = linspace(1,0,length(angPro(t>=116&t<191)));
            angVel = angVMag * angPro;
            t = (0:length(angVel)-1)*dt;
            % Angular acceleration
            angAcc = [zeros(3,1), diff(angVel,[],2)]/dt;
            t = (0:length(angAcc)-1)*dt;
        %% Constant Angular Velocity ON-OFF
        case 2
            %angVMag = [0; 90; 0];
            angVMag = [0; 45; 0];
            Tmo = 160;
            t =-20:dt:Tmo;
            %stop_t = 46.5;
            %stop_t = 60.5;
            %stop_t = 30.5;
            %stop_t = 20.5;
            stop_t = dur;
            angPro = zeros(1,length(t));
            angPro(t>0&t<=1) = linspace(0,1,length(angPro(t>0&t<=1)));
            angPro(t>1&t<stop_t)= 1;
            angPro(t>=stop_t&t<stop_t+1) = linspace(1,0,length(angPro(t>=46.5&t<47.5)));
            angVel = angVMag * angPro;

            % Angular acceleration
            angAcc = [zeros(3,1), diff(angVel,[],2)]/dt;
        %% Angular Velocity Positive-Negative
        case 3
            angVMag = [0; 45; 0];
            Tmo = 90;
            t =-50:dt:Tmo;
            angPro = zeros(1,length(t));
            talpha = 25;
            angPro(t>0.0&t<talpha) = 0.5-0.5*cos(t(t>0.0&t<talpha)*pi/talpha);
            angPro(t>=talpha&t<30)= 1;
            angPro(t>=30&t<31) = 0.5 + 0.5*cos((t(t>=30.0&t<31.0)-30.0)*pi);
            angVel = angVMag * angPro;

            % Angular acceleration
            angAcc = [zeros(3,1), diff(angVel,[],2)]/dt;
        %% Constant Angular Acceleration-Illusion
        case 4
            %angVMag = [0; 360; 0];
            angVMag = [0; 900; 0];
            Tmo = 201;
            t =-50:dt:Tmo;
            angPro = zeros(1,length(t));
            angPro(t>1&t<=76) = linspace(0,1,length(angPro(t>1&t<=76)));
            angPro(t>=76&t<=201)= linspace(1,1,length(angPro(t>=76&t<=201)));
            %angPro(t>=116&t<191) = linspace(1,0,length(angPro(t>=116&t<191)));
            angVel = angVMag * angPro;
            t = (0:length(angVel)-1)*dt;
            % Angular acceleration
            angAcc = [zeros(3,1), diff(angVel,[],2)]/dt;
            t = (0:length(angAcc)-1)*dt;
        %% Karmali Motion Profile
        case 5
            angVMag = [0; 90; 0];
            %angVMag = [0; 45; 0];
            Tmo = 50;
            t =-10:dt:Tmo;
            %stop_t = 46.5;
            stop_t = 61.5;
            %stop_t = 31.5;
            angPro = zeros(1,length(t));
            angPro(t>0.0&t<=1.0) = linspace(0,1,length(angPro(t>0.0&t<=1.0)));
            angPro(t>1.0)= 1;
            angVel = angVMag * angPro;

            % Angular acceleration
            angAcc = [zeros(3,1), diff(angVel,[],2)]/dt;
    end
end