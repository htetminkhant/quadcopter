%Compute torques, given current inputs, length, drag coefficient, and thrust coefficient
function tau = torques(inputs, L, b, k)
    inputs(1) = 2.6e+5;
    inputs(2) = 2.6e+5;
    inputs(3) = 2.6e+5;
    inputs(4) = 2.6e+5;
  tau = [
       L * k * (inputs(1) - inputs(3))
       L * k * (inputs(2) - inputs(4))
       b * (inputs(1) - inputs(2) + inputs(3) - inputs(4))
       ];
end