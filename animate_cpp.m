% if t = 0
t = 0;
dt = 0.05;

%Create plot
createPlot()
%Build to quadcopter.  
[model, thrusts] = quadcopter();
setVisLimits()

%if t > 0
% Simulation times, in seconds.
t = t + dt;
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

function createPlot()
    F = figure(1);
    figure(F)
    plots = [subplot(3, 2, 1:4), subplot(3, 2, 5), subplot(3, 2, 6)];
    subplot(plots(1));
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

function setVisLimits()
    % Set axis scale and labels.
    axis([-10 30 -20 20 5 30]);
    zlabel('Height');
    title('Quadcopter Flight Simulation'); 
end
 
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