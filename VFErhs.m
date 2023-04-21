function rhs = VFErhs(s, XTnb, c0, t)
rhs = zeros (12,1);
%X = XTnb(1:3);
T = XTnb(4:6);
n = XTnb(7:9);
b = XTnb(10:12);

rhs(1:3) = T;
rhs(4:6) = (c0 / sqrt(t)) * n;
rhs(7:9) = - (c0 / sqrt(t)) * T + (s / (2 * t)) * b;
rhs(10:12) = - (s / (2 * t)) * n;
