function rhs = xtyta(t, xy, a)
x = xy(1); %rhs
y = xy(2); %lhs
rhs = zeros(2,1);
rhs(1) = a * y + t;
rhs(2) = a * x - t;
end