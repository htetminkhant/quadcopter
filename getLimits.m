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