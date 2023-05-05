% approximate the derivatives of f(x) = exp(sin(x))
clear
N = 16*1024;
s = 2 * pi * (0:N-1) / N;
f = exp(sin(s));
f_ = fft(f);

f_(abs(f_) / N < eps) = 0; %Krasny's filter

k = [0:N/2-1 -N/2:-1];
fs = ifft(1i * k .* f_);
fss = ifft(1i * k .* 2 .* f_);
fsss = ifft(1i * k .* 3 .* f_);
fssss = ifft(1i * k .* 4 .* f_);
fs10 = ifft(1i * k .* 10 .* f_);
norm(fs - cos(s) .* exp(sin(s)), "inf")
norm(fss - (-sin(s) .* exp(sin(s)) + cos(s) .^ 2 .* exp(sin(s))), "inf")

%vectorize(simplify(diff(exp(sin(s)),3)))
norm(fss - (-(sin(2.*s) .* exp(sin(s)) .* (sin(s) +3))./2), "inf")

%norm(fss - double(subs(simplify(diff(str2sym('exp(sin(s))'),4)),s))

norm(fs10 - double(subs(simplify(diff(str2sym('exp(sin(s))'),10)),s)), "inf")