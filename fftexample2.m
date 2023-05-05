clear
a = 3; b = 2;
N = 1024;
s = 2 * pi * (0:N-1) / N;
x = a * cos(s);
y = b * sin(s);

k = [0:N/2-1 -N/2:-1];
x_ = fft(x);
x_(abs(x_) / N < eps) = 0; %Krasny's filter
xs = real(ifft(1i * k .* x_));

y_ = fft(y);
y_(abs(y_) / N < eps) = 0; %Krasny's filter
ys = real(ifft(1i * k .* y_));

L = (2 * pi / N) * sum(sqrt(xs .^ 2 + ys .^ 2)) %trapezoidal rule

%e = L - 2*pi*3


%plot(x,y);
%axis equal