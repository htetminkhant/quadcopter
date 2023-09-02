%Plot three components of a vector in RGB.
function multiplot(value,time)
    persistent vals times
    % Select the parts of the data to plot.
    times = [times; time];
    vals = [vals; value];
%     if number of elements in vals (length) < 2
%         return
%     end
    % Plot in RGB, with different markers for different components.
    plot(times, value(1, :), 'r-', times, value(2, :), 'g.', times, value(3, :), 'b-.');
    
    % Set axes to remain constant throughout plotting.
    xmin = min(times);
    xmax = max(times) + 1;
    ymin = 1.1 * min(min(vals));
    ymax = 1.1 * max(max(vals));
    axis([xmin xmax ymin ymax]);
end
