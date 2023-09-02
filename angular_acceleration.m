
%Compute angular acceleration in body frame
function omegad = angular_acceleration(inputs, omega, I, L, b, k)
    tau = torques(inputs, L, b, k);
    omegad = inv(I) * (tau - cross(omega, I * omega));
end
