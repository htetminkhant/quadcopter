
 % Visualize plot for motion of quad over time
createplot();


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




 
 

 
