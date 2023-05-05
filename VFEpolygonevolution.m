clear 
close all
M = 3; % number of sides
N = 64 * M; %total number of points
s = 2 * pi * (0:N-1) / N;
T0 = zeros(N, 3);
for m = 0:M-1
    T0((1:N/M) + m * (N/M), 1) = cos(2*pi*m /M);
    T0((1:N/M) + m * (N/M), 2) = sin(2*pi*m /M);
end
X0 = (2 * pi /N) * cumsum([zeros(1,3);T0(1:N-1,:)]);
X0(:,1) = X0(:,1) - mean(X0(:,1));
X0(:,2) = X0(:,2) - mean(X0(:,2));

k = [0:N/2-1 -N/2:-1; 0:N/2-1 -N/2:-1; 0:N/2-1 -N/2:-1].';

X = X0;
T = T0;

%Runge Kutta methods
nmax = 2400;
xcorner = zeros(nmax+1,3);
xcorner(1,:) = X(1,:); 

dt = (2*pi/M^2)/nmax;

XX1 = zeros(N, nmax+1);
XX2 = zeros(N, nmax+1);
XX3 = zeros(N, nmax+1);

TT1 = zeros(N, nmax+1);
TT2 = zeros(N, nmax+1);
TT3 = zeros(N, nmax+1);

XX1(:, m+1) = X(:,1);
XX2(:,m+1) = X(:,2);
XX3(:,m+1) = X(:,3);

TT1(:,1) = T(:,1);
TT1(:,1) = T(:,2);
TT1(:,1) = T(:,3);
for m = 1:nmax
    T_ = fft(T);
    T_(abs(T_) / N < eps) = 0; %Krasny's filter
    Ts = real(ifft(1i * k .* T_));
    Tss = real(ifft((1i * k) .^ 2 .* T_));
    AX = cross(T,Ts);
    AT = cross(T,Tss);
    Taux = T + (dt/2) * AT;

    T_ = fft(Taux);
    T_(abs(T_) / N < eps) = 0; %Krasny's filter
    Ts = real(ifft(1i * k .* T_));
    Tss = real(ifft((1i * k) .^ 2 .* T_));
    BX = cross(Taux,Ts);
    BT = cross(Taux,Tss);
    Taux = T + (dt/2) * BT;

    T_ = fft(Taux);
    T_(abs(T_) / N < eps) = 0; %Krasny's filter
    Ts = real(ifft(1i * k .* T_));
    Tss = real(ifft((1i * k) .^ 2 .* T_));
    CX = cross(T,Ts);
    CT = cross(T,Tss);
    Taux = T + dt * CT;

    T_ = fft(Taux);
    T_(abs(T_) / N < eps) = 0; %Krasny's filter
    Ts = real(ifft(1i * k .* T_));
    Tss = real(ifft((1i * k) .^ 2 .* T_));
    DX = cross(Taux,Ts);
    DT = cross(Taux,Tss);

    X = X + (dt/6) * (AX + 2 * BX + 2 * CX + DX);
    T = T + (dt/6) * (AT + 2 * BT + 2 * CT + DT);

    normT = sqrt(T(:,1) .^ 2 + T(:,2) .^ 2 + T(:,3) .^ 2);
    T = T ./ [normT normT normT];

    %plot3(X(:,1),X(:,2),X(:,3));
    %axis([-2,2,-2,2,-1,1])
    %drawnow

    xcorner(m+1,:) = X(1,:);

    XX1(:, m+1) = X(:,1);
    XX2(:,m+1) = X(:,2);
    XX3(:,m+1) = X(:,3);

    TT1(:,m+1) = T(:,1);
    TT1(:,m+1) = T(:,2);
    TT1(:,m+1) = T(:,3);

    z = -sqrt(xcorner(:,1).^2 + xcorner(:,2).^2 + 1i*xcorner(:,3).^2);
    plot(z(1:m+1))
    drawnow

end

%plot3(T(:,1),T(:,2),T(:,3), '.');
%plot3(xcorner(:,1),xcorner(:,2),xcorner(:,3))
%plot3(X(:,1),X(:,2),X(:,3), '.')
%plot(X(:,1),X(:,2)), axis equal