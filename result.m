% Initial system state.
x = [0; 0; 10];
theta = zeros(3,1);
xdot = zeros(3,1);

% Deviation is in degrees/sec.
deviation = 100;
thetadot = deg2rad(2 * deviation * rand(3,1) - deviation); 

F = figure(1);
figure(F)
plots = [subplot(3, 2, 1:4), subplot(3, 2, 5), subplot(3, 2, 6)];
subplot(plots(1));

%Build to quadcopter.  
[model, thrusts] = quadcopter();

% Set axis scale and labels.
axis([-10 30 -20 20 5 30]);
zlabel('Height');
title('Quadcopter Flight Simulation');    

% Animate the quadcopter with data from the simulation.
animate(x,theta,xdot,thetadot,model,thrusts,plots) 
