% Physical constants.
g = 9.81;
m = 0.5;
L = 0.25;
k = 3e-6;
b = 1e-7;
I = diag([5e-3, 5e-3, 10e-3]);
kd = 0.25;
% Simulation times, in seconds.
Start_time = 0;
End_time = 10;
dt = 0.05;
times = Start_time:dt:End_time;
% Number of points in the simulation.
N = numel(times);
% Initial system state.
x = [0; 0; 10];
xdot = zeros(3,1);
theta = zeros(3,1);
% Deviation is in degrees/sec.
deviation = 100;
thetadot = deg2rad(2 * deviation * rand(3,1) - deviation);
 
ff = figure(1);
figure(ff)
plots = [subplot(3, 2, 1:4), subplot(3, 2, 5), subplot(3, 2, 6)];

for t =  times 
    subplot(plots(1));
    %Build to quadcopter.  
    [model, thrusts] = quadcopter();
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
    % Compute translation to correct linear position coordinates.
    move = makehgtform('translate', x);
    % Compute rotation to correct angles. Then, turn this rotation into a 4x4 matrix represting this affine transformation.
    angles = theta;
    rotate = rotation(angles);
    rotate = [rotate, zeros(3, 1); zeros(1, 3), 1];
    % Move the quadcopter to the right place, after putting it in the correct orientation.
    set(model,'Matrix', move * rotate);
    % Compute scaling for the thrust cylinders.
    scales = exp(inputvalue()/min(abs(inputvalue())) + 5) - exp(6) + 1.5;
    for i = 1:4
        s = scales(i);
        if s < 0
            scalez = makehgtform('yrotate', pi)  * makehgtform('scale', [1, 1, abs(s)]);
        elseif s > 0
            scalez = makehgtform('scale', [1, 1, s]);
        end
        set(thrusts(i), 'Matrix', move * rotate * scalez);
    end
    % Set axis scale and labels.
    axis([-10 30 -20 20 5 30]);
    zlabel('Height');
    title('Quadcopter Flight Simulation'); 
    % Update the drawing.        
    xmin = x(1)-20;
    xmax = x(1)+20;
    ymin = x(2)-20;
    ymax = x(2)+20;
    zmin = x(3)-5;
    zmax = x(3)+5;
    axis([xmin xmax ymin ymax zmin zmax]);
    hold off;
    drawnow;
    % Animate the quadcopter with data from the simulation.
    animate(x,theta,xdot,thetadot,model,thrusts);
    % Use the bottom two parts for angular velocity and displacement.   
    subplot(plots(2));
    multiplot(thetadot,t);
    xlabel('Time (s)');
    ylabel('Angular Velocity (rad/s)');
    title('Angular Velocity');

    subplot(plots(3));
    multiplot(theta,t);
    xlabel('Time (s)');
    ylabel('Angular Displacement (rad)');
    title('Angular Displacement');
 end


