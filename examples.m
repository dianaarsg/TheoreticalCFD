syms x y a z(t)
format long
int(1/(1+x^2)^10) %Integration
diff(sin(2*x+3*y+4*a),x)
diff(sin(2*x+3*y+4*a),y)
diff(sin(2*x+3*y+4*a),a)
diff(sin(2*x+3*y+4*a),x,100)
diff(sin(2*x+3*y+4*a),x,x,x,y,y,y,a,a)
int(exp(-x^2), x, -inf, inf)
int(x^10*exp(-x^2), x, -inf, inf)

solve(x+y == 3, x-y ==1) %sol is a strap to retrieve the solution

[solx, soly] = solve (x + y == 3 + a, x - y == 1 - a, x, y)

sol = solve(x^2 + y^2 == 3, x^2 - y^2 == 1)

double(solve(x^3 + x + 1))

%ode45 integrator
dsolve(diff(z(x))== z(x))

dsolve(diff(z)== z)

dsolve(diff(z,2)== z, z(0) ==1, subs(diff(z), 0) ==3)
dsolve(diff(z)==z^2, z(0) == 1)


%vpa
vpa(pi, 1000)
%str2sym 
vpa(x + str2sym('pi^2'), 100)

vpa(double(solve(x^3 + x + 1)),100)