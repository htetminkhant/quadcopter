% Plot three components of a vector in RGB.
function multiplot(time,value)
    persistent vals times
    
    % Select the parts of the data to plot.
    
    times = [times; time];
    vals =  [vals; value];
     
    % Plot in RGB, with different markers for different components.
    plot(times,vals(1),'r-',times,vals(2),'g.',times,vals(3),'b.' );
    hold on;
    
    % Set axes to remain constant throughout plotting.
    xmin = min(times);
    xmax = max(times)+1;
    ymin = 1.1 * min(min(vals));
    ymax = 1.1 * max(max(vals));
    axis([xmin xmax ymin ymax]);
end
