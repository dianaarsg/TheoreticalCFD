clear 
syms t x(t) y(t)
sol = dsolve(diff(x) == y + t, diff(y) == x - t, x(0) == 2, y(0) == 3);
xsol = simplify(sol.x)
ysol = simplify(sol.y)
options  = odeset ('RelTol',3e-14, 'AbsTol',1e-15);
[T,XY] = ode45(@xtyt, 0:0.01:2, [2; 3], options);
X = XY(:,1);
Y = XY(:,2);
norm(X - double(subs(xsol, T)), "inf")
norm(Y - double(subs(ysol, T)), "inf")

%% VFE
clear

v = VideoWriter("VFEmovie.mp4", "MPEG-4");
v.Quality = 100; %fps
open(v);

c0 = 1;
t = 5;
L = 50;
T0 = [1 0 0];
n0 = [0 1 0];
b0 = [0 0 1];


options  = odeset ('RelTol',3e-14, 'AbsTol',1e-15);
for t = 1e-2:1e-2:10
    X0 = 2 * c0 * sqrt(t) * [0, 0, 1];
    XTnb0 = [X0 T0 n0 b0];
    VFErhsaux = @(s , XTnb)VFErhs(s, XTnb, c0, t);
    [spos, XTnbpos] = ode45(VFErhsaux, 0:0.01:L, XTnb0);
    [sneg, XTnbneg] = ode45(VFErhsaux, 0:-0.01:-L, XTnb0);
    S = [sneg(end: -1 : 2); spos];
    XTnb = [XTnbneg(end:-1:2, :);XTnbpos];
    X = XTnb(:, 1:3);
    T = XTnb(:, 4:6);
    plot3(X(:,1), X(:,2), X(:,3));
    %plot3(T(:,1),T(:,2), T(:,3));
    title(['t = ', num2str(t)]);
    writeVideo(v, getframe(gcf));
    drawnow %to plot until is done
end

close (v)
%figure;
%plot3(T(:,1),T(:,2), T(:,3));
%drawnow;

%% Finite differences
clear
close all
tic

c0 = 1;
t0 = 1;
tend = 0;
nmax = 10000;
dt = (t0 - tend)/nmax;
L = 20;
ds = 0.01/10;
X0 = 2 * c0 * sqrt(t0) * [0, 0, 1];
T0 = [1 0 0];
n0 = [0 1 0];
b0 = [0 0 1];
XTnb0 = [X0 T0 n0 b0];

options  = odeset ('RelTol',3e-14, 'AbsTol',1e-15);
VFErhsaux = @(s , XTnb)VFErhs(s, XTnb, c0, t0);
[spos, XTnbpos] = ode45(VFErhsaux, 0:ds:L, XTnb0, options);
[sneg, XTnbneg] = ode45(VFErhsaux, 0:-ds:-L, XTnb0, options);

S = [sneg(end: -1 : 2); spos];
XTnb = [XTnbneg(end:-1:2, :);XTnbpos];

%Initial data
XX0 = XTnb(:, 1:3);
TT0 = XTnb(:, 4:6);

XX = XX0;
TT = TT0;
%XXs = (XX(3:end, :) - XX(1:end-2, :)) / (2 * ds);
%max(max(abs(XXs - TT(2:end-1,:))));

%Runge-Kutta method
% h = dt
for m = 1:nmax
    XXs = [zeros(1,3); (XX(3:end, :) - XX(1:end-2, :)) / (2 * ds);zeros(1,3)]; %first derivative 
    TTss = [zeros(1,3);(XX(3:end, :) - 2 * XX(2:end-1, :) + XX(1:end-2, :)) / ds^2; zeros(1,3)]; %2nd derivative
    k1X = cross(XX, XXs);
    %k1T = cross(XX, XXss);

    %XXaux = XX + 0.5 * dt * k1X;
    XXaux = XX + 0.5 * dt * k1X;
    XXs = [zeros(1,3); (XXaux(3:end, :) - XXaux(1:end-2, :)) / (2 * ds);zeros(1,3)]; %first derivative 
    XXss = [zeros(1,3);(XXaux(3:end, :) - 2 * XXaux(2:end-1, :) + XXaux(1:end-2, :)) / ds^2; zeros(1,3)]; %2nd derivative
    k2X = cross(XXaux, XXs);
    %k2X = cross(XXaux, XXss);

    XXaux = XX + 0.5 * dt * k2X;
    XXs = [zeros(1,3); (XXaux(3:end, :) - XXaux(1:end-2, :)) / (2 * ds); zeros(1,3)]; %first derivative 
    TTss = [zeros(1,3);(XXaux(3:end, :) - 2 * XXaux(2:end-1, :) + XXaux(1:end-2, :)) / ds^2; zeros(1,3)]; %2nd derivative
    k3X = cross(XXaux, XXs);
    %k3T = cross(TTaux, TTss);

    XXaux = XX + dt * k3X;
    XXs = [zeros(1,3); (XXaux(3:end, :) - XXaux(1:end-2, :)) / (2 * ds);zeros(1,3)]; %first derivative 
    XXss = [zeros(1,3);(XXaux(3:end, :) - 2 * XXaux(2:end-1, :) + XXaux(1:end-2, :)) / ds^2; zeros(1,3)]; %2nd derivative
    k4X = cross(XXaux, XXss);
    %k4T = cross(TTaux, TTss);

    XX = XX + (dt/6) * (k1X + 2 * k2X + 2 * k3X + k4X);
    %TT = TT + (dt/6) * (k1T + 2 * k2T + 2 * k3T + k4T);
    if isnan(XX(2,:)), error ('Stability error!'); end

    plot3(X(:,1), X(:,2), X(:,3));
    %plot3(T(:,1),T(:,2), T(:,3));
    %title(['t = ', num2str(t)]);
    drawnow %to plot until is done

end

toc

%plot3(X(:,1), X(:,2), X(:,3));
%plot3(T(:,1),T(:,2), T(:,3));
%title(['t = ', num2str(t)]);




