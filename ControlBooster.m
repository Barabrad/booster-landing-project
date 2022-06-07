function [angAtt, phase, mdot, Fx] = ControlBooster(t, DSuicide, state, parameters)
    global p4Y time2 time3 oldMdot oldFx;
% ASTE 101, Fall 2020

    g = 9.8;

    h = parameters(1);
    D = parameters(2);
    Ms = parameters(3);
    Isp = parameters(4);
    mdotMax = parameters(5);
    FxMax = parameters(6);

    x = state(1);
    y = state(2);
    Vx = state(3);
    Vy = state(4);
    theta = state(5);
    omega = state(6);
    M = state(7);
    I = rocketMoment(M);
    
    ueq = Isp * g;
    
    %Take cos then acos to "synchronise" angles
    angleR = cosd(theta - 180); %Angle of the rocket relative to vertical
    angleP = cosd(90 - atan2d(Vy, Vx)); %Angle of path rel. to ver.
    angAtt = acosd(angleR) - acosd(angleP);
    
    if abs(Vx) < 10
        thetaLQR = theta; %Goal: zero out theta
    else
        thetaLQR = angAtt; %Zeroing out angle of attack reduces Vx
    end
    
    function I = rocketMoment( M )
        I = (M/12) * ( 3*(D/2)^2 + h^2 );
    end

    function Fx = CorrectAngle(angleDiff, omega)
        LQRs3 = [ angleDiff, omega ]';

        A = [ 0 1;...
            0 0 ];

        B = [ 0 (0.5*h/I)*(180/pi) ]';
        
        Q = diag( [ 2 5 ] );

        R = 1e-6;
        N = 0;

        k = lqr( A, B, Q, R, N );
        u = -k * LQRs3;

        Fx = u;
    end

% Phase 1: End-over-end
    tEOEhalf = ceil(sqrt( 2*pi*I/(FxMax*h) ));
    if t < tEOEhalf
        phase = "End-over-end (a)";
        mdot = 0;
        Fx = FxMax;
        return
    elseif (t >= tEOEhalf && t < 2*tEOEhalf)
        phase = "End-over-end (d)";
        mdot = 0;
        Fx = -FxMax;
        return
    elseif (t >= 2*tEOEhalf && t < 2*tEOEhalf + 0.2)
        phase = "Stopping Fx";
        mdot = 0;
        Fx = 0;
        return
    end
    
% Phase 2: Coast
    if (y > DSuicide)
        phase = "Coasting";
        mdot = -x/6300; %A small arbitrary force to reduce overshot of x=0
        Fx = CorrectAngle(thetaLQR, omega);
        time2 = t;
        time3 = t;
    else
% Phase 3: Decelerate
        landingVy = -100;
        if (Vy <= landingVy && time2 == time3)
            phase = "Slowing";
            mdot = mdotMax;
            Fx = CorrectAngle(thetaLQR, omega);
            time2 = t;
            time3 = t;
        elseif (Vy > -1)
            phase = "Landing";
            mdot = 0;
            Fx = CorrectAngle(thetaLQR, omega);
        else
% Phase 4: Landing
            if time2 == time3
                p4Y = y; %The height of switch from dec. to landing
            end
            phase = "Landing";
            %We need a transition from Vy0 to Vy = -1
            buffer = 30;
            if y > buffer
                VDesired = landingVy + (1-landingVy)*(p4Y-y)/(p4Y-buffer);
            else
                VDesired = -1;
            end
            
            LQRs2 = [ x, Vx, thetaLQR, omega, y, Vy-VDesired ]';
            
            if oldMdot <= 0
                Fth0 = M*g;
            else
                Fth0 = oldMdot*ueq;
            end

            A = [ 0 1 0 0 0 0;...
                0 0 (Fth0/M)*(pi/180) 0 0 0;...
                0 0 0 1 0 0;...
                0 0 0 0 0 0;...
                0 0 0 0 0 1;...
                0 0 (oldFx/M)*(pi/180) 0 0 0];

            B = [ 0 (1/M) 0 (0.5*h/I)*(180/pi) 0 0;...
                0 0 0 0 0 (ueq/M)]';

            Q = diag( [ 1 2 5 10 1 5 ] );

            R = [1e-6 0;...
                0 1e-6];
            N = [0 0; 0 0; 0 0; 0 0; 0 0; 0 0];

            k = lqr( A, B, Q, R, N );
            u = -k * LQRs2;

            Fx = u(1);
            
            if u(2) < 0
                mdot = M/Isp; %keeps Vy constant
            else
                mdot = u(2) + M/Isp;
            end
            
            if (Vy < landingVy-50) %Switch back to dec. if speed increases
                time2 = t;
                time3 = t;
            else
                time3 = t;
            end
        end
    end
    oldMdot = mdot;
    oldFx = Fx;
end
