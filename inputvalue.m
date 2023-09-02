function in = inputvalue()
    in = zeros(4, 1);
    in(:) = 550;
    in(1) = in(1) + 168;
    in(3) = in(3) + 168;
    in = in .^ 2;
end