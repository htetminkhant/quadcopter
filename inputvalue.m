function in = inputvalue()
    in = zeros(4, 1);
    in(:) = 7.5e+5;
    in(1) = in(1) + 150;
    in(3) = in(3) + 150;
    in = in .^ 2;
end