function simulateBooster(controlBooster)
% ASTE 101, Fall 2020
% Main program for running booster landing simulation

    if nargin < 1
        controlBooster = @ControlBooster;
    end

    % Physical quantities
    g = 9.8; % Earth gravity (constant, m/s^2)
    km = 1000; 

    % Assumed quantities
    h = 50; % Height of booster, m
    D = 5; % Diameter of booster, m
    Ms = 1e4; % Empty (structural) mass of booster, kg
    R = 12; % Mass ratio after BECO
    Isp = 350; % Specific impulse of rocket engines
    alphaMax0 = 5; % Max. angular acceleration (deg/s^2) after BECO
    accelMax0 = 2*g; % Acceleration at maximum mdot after BECO
    mdotError = 0.01; % Commanded mass flow has this much fractional error
    FxError = 0.05; % Commanded sideways force has this much fractional error

    % For computation
    dt = 0.1; % Timestep

    % Derived quantities
    Mp0 = (R - 1) * Ms; % Propellant mass after BECO
    M0 = Ms + Mp0; % Rocket mass after BECO
    mdotMax = accelMax0 / (Isp*g) * R*Ms;
    I0 = (M0/12) * ( 3*(D/2)^2 + h^2 ); % Initial moment of inertia
    % Fx h/2 = I0 * alpha
    FxMax = I0 * alphaMax0 / (h/2) * pi/180;
    ueq = Isp*g;

    % Initial conditions
    % Per https://www.reddit.com/r/spacex/comments/bcdss4/some_meco_and_beco_speeds/
    % BECO cuts off 1611 m/s, altitude 58 km, T=02:35
    x0 = -200 * km;
    y0 = 60 * km;
    V0 = 1600;
    theta0 = 80;
    omega0 = 0;
    Vx0 = V0 * sind(theta0);
    Vy0 = V0 * cosd(theta0);

    state0 = [ x0 y0 Vx0 Vy0 theta0 omega0 M0 ]';
    state = state0;

    parameters = [ h, D, Ms, Isp, mdotMax, FxMax ];

    t = 0;

    function I = rocketMoment(M)
        I = (M/12) * ( 3*(D/2)^2 + h^2 );
    end

    function q = addError(qNominal, fractionalError)
        mu = qNominal;
        sigma = abs(fractionalError * mu);
        q = normrnd(mu, sigma);
    end

    %Found this on Desmos, but I added the buffer multiplier
    function D = suicideBurn(Vy, M, g, y, ueq, mdot, bufferMult)
        T = M*(1-exp( (Vy - sqrt(2*g*y))/(ueq) ))/mdot;
        D = (T/2)*(-Vy + sqrt(2*g*y))*bufferMult;
    end

    close all;
    figure;

    %Height at which to start suicide burn
    DSuicide = suicideBurn(Vy0, M0, g, y0, ueq, mdotMax, 3);

    % Main simulation loop
    i = 1;
    fprintf( 'Initial: t=%f, x=%f, y=%f, Vx=%f, Vy=%f, theta=%f, omega=%f\n', ...
        t, x0, y0, Vx0, Vy0, theta0, omega0 );
    while( state(2) > 0 && t < 1000 )
        [angAtt, phase, mdot, Fx] = controlBooster(t, DSuicide, state, parameters);
        Fx = max( -FxMax, min( FxMax, Fx ) ); % Ensure legal value
        Fx = addError( Fx, FxError ); % Add noise
        x = state(1);
        y = state(2);
        Vx = state(3);
        Vy = state(4);
        theta = -180 + mod( state(5)+180, 360 );
        omega = state(6);
        M = state(7);
        I = rocketMoment(M);
        if( M > Ms )
            mdot = max( 0, min( mdotMax, mdot ) ); % Ensure legal value
            mdot = addError(mdot, mdotError); % Add noise
            accel = mdot*ueq / M;
            ax = accel*sind(theta) + Fx*cosd(theta)/M;
            ay = accel*cosd(theta) - Fx*sind(theta)/M;
            Vx = Vx + ax*dt;
            Vy = Vy + ay*dt;
            M = M - mdot*dt;
            if M <= Ms
                fprintf("/nFuel spent!/n");
            end
        else
            mdot = 0;
            Fx = 0;
        end
        Vy = Vy - g*dt; % Include gravity
        x = x + Vx*dt;
        y = y + Vy*dt;

        alpha = Fx * h/2 / I * 180/pi;
        theta = theta + omega*dt;
        omega = omega + alpha*dt;

        state = [ x, y, Vx, Vy, theta, omega, M ]';

        tvals(i) = t;
        xvals(i) = x;
        yvals(i) = y;
        Vxvals(i) = Vx;
        Vyvals(i) = Vy;
        thetavals(i) = theta;
        omegavals(i) = omega;
        attackvals(i) = angAtt;
        Mvals(i) = M;

        t = t + dt;
        i = i + 1;

        fprintf( '%s: mdot=%f, Fx=%f, t=%f, x=%f, y=%f, Vx=%f, Vy=%f, theta=%f, omega=%f\n', ...
            phase, mdot, Fx, t, x, y, Vx, Vy, theta, omega );

    end
    x = state(1);
    y = state(2);
    Vx = state(3);
    Vy = state(4);
    theta = state(5);
    omega = state(6);
    
    if (Vy < -5 || abs(Vx) > 2 || abs(theta) > 5 || abs(omega) > 10)
        keyTerm = "crashed";
    else
        keyTerm = "landed";
    end

    fprintf( '\nBooster %s at t=%f, x=%f, Vx=%f, Vy=%f, theta=%f, omega=%f\n', ...
        keyTerm, t, x, Vx, Vy, theta, omega );
    fprintf( 'Propellant Used: %f kg out of %f kg\n', M0 - M, Mp0 );

    subplot(4,2,1);
    plot(tvals, xvals);
    ylabel( 'x, m' );

    subplot(4,2,2);
    plot(tvals, yvals);
    ylabel( 'y, m' );

    subplot(4,2,3);
    plot(tvals, Vxvals);
    ylabel( 'Vx, m/s' );

    subplot(4,2,4);
    plot(tvals, Vyvals);
    ylabel( 'Vy, m/s' );

    subplot(4,2,5);
    plot(tvals, thetavals);
    ylabel( '\theta, degrees' );

    subplot(4,2,6);
    plot(tvals, attackvals);
    ylabel( 'Angle of attack, degrees' );
    
    subplot(4,2,7);
    plot(tvals, omegavals);
    ylabel( '\omega, degrees/s' );

    subplot(4,2,8);
    plot(tvals, Mvals);
    ylabel( 'Rocket mass M, kg' );
    xlabel( 'Time, s' );

end
  
  
  

