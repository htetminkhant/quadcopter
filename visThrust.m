function visThrust(thrusts, move, rotate, scale)
     % Compute scaling for the thrust cylinders.
    scales = exp(scale/min(abs(scale)) + 5) - exp(6) + 1.5;
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