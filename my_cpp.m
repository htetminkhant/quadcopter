function my_cpp(x,theta,xdot,thetadot,model,thrusts,plots)

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
visVelocities(t, theta, thetadot)
end
end 


function [rotate] = getRotation(theta)
    % Compute rotation to correct angles. Then, turn this rotation into a 4x4 matrix represting this affine transformation.
    rotate = rotation(theta);
    rotate = [rotate, zeros(3, 1); zeros(1, 3), 1];
end
% end function getRotation

 function [xmin, xmax, ymin, ymax, zmin, zmax] = getLimits(x)
    % Update the drawing.      
    xmin = x(1)-20;
    xmax = x(1)+20;
    ymin = x(2)-20;
    ymax = x(2)+20;
    zmin = x(3)-5;
    zmax = x(3)+5;
 end
 % end function getLimits
 
function visThrust(thrusts, move, rotate)
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
 end
 % end function visThrust
 
 function visVelocities(t, theta, thetadot)
     % Use the bottom two parts for angular velocity and displacement.   
    subplot(plots(2));
    multiplot(t,thetadot);
    xlabel('Time (s)');
    ylabel('Angular Velocity (rad/s)');
    title('Angular Velocity');

    subplot(plots(3));
    multiplot(t,theta);
    xlabel('Time (s)');
    ylabel('Angular Displacement (rad)');
    title('Angular Displacement'); 
 end