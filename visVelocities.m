 function visVelocities(t, theta, thetadot)
     % Use the bottom two parts for angular velocity and displacement.
    plots = createplot();
    subplot(plots(2));
    multiplot(t,thetadot);
    xlabel('Time (s)');
    ylabel('Angular Velocity (rad/s)');
    title('Angular Velocity');
    
    % Use the bottom two parts for angular velocity and displacement.
    subplot(plots(3));
    multiplot(t,theta);
    xlabel('Time (s)');
    ylabel('Angular Displacement (rad)');
    title('Angular Displacement');

 end