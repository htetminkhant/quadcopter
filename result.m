% Initial system state.
x = [0; 0; 10];
theta = zeros(3,1);
xdot = zeros(3,1);

% Deviation is in degrees/sec.
deviation = 100;
thetadot = deg2rad(2 * deviation * rand(3,1) - deviation); 

createfigure();

plots = createplot();
subplot(plots(1));


%Build to quadcopter.  
[model, thrusts] = quadcopter();

createlabel();   

% Animate the quadcopter with data from the simulation.
% Simulation times, in seconds.
Start_time = 0;
End_time = 4;
dt = 0.05;
times = Start_time:dt:End_time;

ind = 0;
for t =  times   
    ind = ind + 1;   
    % Physical constants.
    g = 9.81;
    m = 0.5;
    L = 0.25;
    k = 3e-6;
    b = 1e-7;
    I = diag([5e-3, 5e-3, 10e-3]);
    kd = 0.25;
   
    % Compute forces, torques, and accelerations.
    omega = thetadot2omega(thetadot, theta);
    a = acceleration(t, theta, xdot, m, g, k, kd);
    omegadot = angular_acceleration(t, omega, I, L, b, k);
   
    % Advance system state.
    omega = omega + dt * omegadot;
    thetadot = omega2thetadot(omega, theta);
    theta = theta + dt * thetadot;
    xdot = xdot + dt * a;
    x = x + dt * xdot;
    
   % Visualize plot for motion of quad over time
    plots = createplot();
    subplot(plots(1));


    % Compute translation to correct linear position coordinates.
    move = makehgtform('translate',x);

    rotate = getRotation(theta);    
    % Move the quadcopter to the right place, after putting it in the correct orientation.
    set(model,'Matrix', move * rotate);

    visThrust(thrusts, move, rotate)

    [xmin, xmax, ymin, ymax, zmin, zmax] = getLimits(x);
    axis([xmin xmax ymin ymax zmin zmax]);
    %hold off;
    drawnow;
    visVelocities(t, theta, thetadot);
end 
