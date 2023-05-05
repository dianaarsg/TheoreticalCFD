clc;
clear all;
[T,X] = ode45(@(t,x)x,[0,10],1);
%x'(t) = f(t, x(t))
norm(X -exp(T), "inf")

options  = odeset ('RelTol',3e-15, 'AbsTol',1e-15);
[T,X] = ode45(@(t,x)x,[0,10],1, options);
norm(X -exp(T), "inf")
%plot(T,X)

options  = odeset ('RelTol',3e-14, 'AbsTol',1e-15);
[X,Y] = ode45(@(t,x)x,[0,10],1, options);
%norm(X -exp(T), "inf")

%syms x y(x), dsolve(diff(y) == y + x, y(0) ==1)

a = 3;
options  = odeset ('RelTol',3e-14, 'AbsTol',1e-15);
[X,Y] = ode45(@(x,y)a*y,[0,1] ,1, options);
norm(Y - exp(a*X), "inf")

clear
syms t x(t) y(t)
sol = dsolve(diff(x) == y + t, diff(y) == x, x(0) == 2, y(0) == 3);
X = simplify(sol.x)
Y = simplify(sol.y)

