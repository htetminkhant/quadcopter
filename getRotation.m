function [rotate] = getRotation(theta)
    % Compute rotation to correct angles. Then, turn this rotation into a 4x4 matrix represting this affine transformation.
    angles = theta;
    rotate = rotation(angles);
    rotate = [rotate, zeros(3, 1); zeros(1, 3), 1];
end
% end function getRotation
