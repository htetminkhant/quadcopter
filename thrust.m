%Compute thrust given current inputs and thrust coefficient.
function T = thrust(inputs, k)
%   inputs(1) = 7.5e+5;
%   inputs(2) = 7.5e+5;
%   inputs(3) = 7.5e+5;
%   inputs(4) = 7.5e+5;
  T = [0; 0; k*sum(inputs)];
end