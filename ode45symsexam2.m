clear 
syms a t x(t) y(t)
sol = dsolve(diff(x) == a * y + t, diff(y) == a * x - t, x(0) == 2, y(0) == 3), t;
xsol = simplify(sol.x)
ysol = simplify(sol.y)

a = 0.5;
anum = 0.5;
options  = odeset ('RelTol',3e-14, 'AbsTol',1e-15);
xtyt = @(t, xy)xtyta(t, xy, anum);
[T,XY] = ode45(xtyt, 0:0.01:2, [2; 3], options);
X = XY(:,1);
Y = XY(:,2);
norm(X - double(subs(subs(xsol, T), a, anum)), "inf")
norm(Y - double(subs(subs(ysol, T), a, anum)), "inf")

%%
