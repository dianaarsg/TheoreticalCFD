% EXERCISE 1 - Runge-Kutta scheme
clear
close all
tic
c0 = 0.2;
t = 1;
L = 20;
ds = 0.05; % ds = 5e-2;
X0 = 2 * c0 * sqrt(t) * [0 0 1];
T0 = [1 0 0];
n0 = [0 1 0];
b0 = [0 0 1];
% whos in command window  to see the defined parameters

mmax = round(L / ds); %rounding result to make it storable by Matlab
X = X0; %define intitial conditions?
T = T0;
n = n0;
b = b0;
s = 0;
XXpos = zeros(mmax+1, 3); XX(1, :) = X0;
TTpos = zeros(mmax+1, 3); TT(1, :) = T0;
nnpos = zeros(mmax+1, 3); nn(1, :) = n0;
bbpos = zeros(mmax+1, 3); bb(1, :) = b0;
%implement runge-kutta method
for m = 1:mmax 
    k1X = T;
    k1T = (c0 / sqrt(t)) * n;
    k1n = -(c0 / sqrt(t)) * T + (s / (2 * t)) * b;
    k1b = -(s / (2 * t)) * n;

    saux = s + .5 * ds;
    %Xaux = X + .5 * ds * k1X; %not used anyway in the right hand side of
    %the function (see notes)
    Taux = T + .5 * ds * k1T;
    naux = n + .5 * ds * k1n;
    baux = b + .5 * ds * k1b;

    k2X = Taux;
    k2T = (c0 / sqrt(t)) * naux;
    k2n = -(c0 / sqrt(t)) * Taux + (saux / (2 * t)) * baux;
    k2b = -(saux / (2 * t)) * naux;

    %saux = s + .5 * ds; %it is the same as above, so it doesnt have to be
    %defined again
    %Xaux = X + .5 * ds * k2X;
    Taux = T + .5 * ds * k2T;
    naux = n + .5 * ds * k2n;
    baux = b + .5 * ds * k2b;

    k3X = Taux;
    k3T = (c0 / sqrt(t)) * naux;
    k3n = -(c0 / sqrt(t)) * Taux + (saux / (2 * t)) * baux;
    k3b = -(saux / (2 * t)) * naux;

    %saux = s + .5 * ds;
    %Xaux = X + .5 * ds * k3X;
    Taux = T + .5 * ds * k3T;
    naux = n + .5 * ds * k3n;
    baux = b + .5 * ds * k3b;

    k4X = Taux;
    k4T = (c0 / sqrt(t)) * naux;
    k4n = -(c0 / sqrt(t)) * Taux + (saux / (2 * t)) * baux;
    k4b = -(saux / (2 * t)) * naux;

    s = m * ds; % s = s + ds;
    X = X + (ds/6) * (k1X + 2 * k2X + 2 * k3X + k4X);
    T = T + (ds/6) * (k1T + 2 * k2T + 2 * k3T + k4T);
    n = n + (ds/6) * (k1n + 2 * k2n + 2 * k3n + k4n);
    b = b + (ds/6) * (k1b + 2 * k2b + 2 * k3b + k4b);

   % T = T / norm(T);
   % n = n / norm(n);
   % b = b / norm(b);

 XXpos(m+1, :) = X;
 TTpos(m+1, :) = T;
 nnpos(m+1, :) = n;
 bbpos(m+1, :) = b;

end

ds = -ds

XXneg = zeros(mmax+1, 3); XX(1, :) = X0;
TTneg = zeros(mmax+1, 3); TT(1, :) = T0;
nnneg = zeros(mmax+1, 3); nn(1, :) = n0;
bbneg = zeros(mmax+1, 3); bb(1, :) = b0;
%implement runge-kutta method
for m = 1:mmax 
    k1X = T;
    k1T = (c0 / sqrt(t)) * n;
    k1n = -(c0 / sqrt(t)) * T + (s / (2 * t)) * b;
    k1b = -(s / (2 * t)) * n;

    saux = s + .5 * ds;
    %Xaux = X + .5 * ds * k1X; %not used anyway in the right hand side of
    %the function (see notes)
    Taux = T + .5 * ds * k1T;
    naux = n + .5 * ds * k1n;
    baux = b + .5 * ds * k1b;

    k2X = Taux;
    k2T = (c0 / sqrt(t)) * naux;
    k2n = -(c0 / sqrt(t)) * Taux + (saux / (2 * t)) * baux;
    k2b = -(saux / (2 * t)) * naux;

    %saux = s + .5 * ds; %it is the same as above, so it doesnt have to be
    %defined again
    %Xaux = X + .5 * ds * k2X;
    Taux = T + .5 * ds * k2T;
    naux = n + .5 * ds * k2n;
    baux = b + .5 * ds * k2b;

    k3X = Taux;
    k3T = (c0 / sqrt(t)) * naux;
    k3n = -(c0 / sqrt(t)) * Taux + (saux / (2 * t)) * baux;
    k3b = -(saux / (2 * t)) * naux;

    %saux = s + .5 * ds;
    %Xaux = X + .5 * ds * k3X;
    Taux = T + .5 * ds * k3T;
    naux = n + .5 * ds * k3n;
    baux = b + .5 * ds * k3b;

    k4X = Taux;
    k4T = (c0 / sqrt(t)) * naux;
    k4n = -(c0 / sqrt(t)) * Taux + (saux / (2 * t)) * baux;
    k4b = -(saux / (2 * t)) * naux;

    s = m * ds; % s = s + ds;
    X = X + (ds/6) * (k1X + 2 * k2X + 2 * k3X + k4X);
    T = T + (ds/6) * (k1T + 2 * k2T + 2 * k3T + k4T);
    n = n + (ds/6) * (k1n + 2 * k2n + 2 * k3n + k4n);
    b = b + (ds/6) * (k1b + 2 * k2b + 2 * k3b + k4b);

   % T = T / norm(T);
   % n = n / norm(n);
   % b = b / norm(b);

 XXneg(m+1, :) = X;
 TTneg(m+1, :) = T;
 nnneg(m+1, :) = n;
 bbneg(m+1, :) = b;

end

toc

XX = [XXneg(end: -1:2, :); XXpos]
TT = [TTneg(end: -1:2, :); TTpos]
nn = [nnneg(end: -1:2, :); nnpos]
bb = [bbneg(end: -1:2, :); bbpos]

% in command line: plot3(XX(:, 1), XX(:, 2), XX(:, 3))
%plot3(XXpos(:, 1), XXpos(:, 2), XXpos(:, 3))
%figure
%plot3(TTpos(:, 1), TTpos(:, 2), TTpos(:, 3))
%figure
%plot3(XXneg(:, 1), XXneg(:, 2), XXneg(:, 3))
%figure
%plot3(TTneg(:, 1), TTneg(:, 2), TTneg(:, 3))
%figure

plot3(XX(:, 1), XX(:, 2), XX(:, 3))
figure
plot3(TT(:, 1), TT(:, 2), TT(:, 3))
figure
plot3(nn(:, 1), nn(:, 2), nn(:, 3))
figure
plot3(bb(:, 1), bb(:, 2), bb(:, 3))
figure