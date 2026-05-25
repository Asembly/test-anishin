pkg load signal;

Fs = 100;
A = 1;
tau = 2;
T = 3;
Tl = max(tau, T) * 2;
t = -Tl:1/Fs:Tl;

s_rect = zeros(size(t));
s_rect(t >= 0 & t <= tau) = A;

s_tri1 = zeros(size(t));
s_tri1(t >= 0 & t <= T) = A * (t(t >= 0 & t <= T) / T);

s_tri2 = zeros(size(t));
s_tri2(abs(t) <= T) = A * (1 - abs(t(abs(t) <= T)) / T);

c7 = xcorr(s_rect, s_tri1);
c7 = c7 / max(abs(c7));

c8 = xcorr(s_rect, s_tri2);
c8 = c8 / max(abs(c8));

t_c = (0:length(c7)-1)/Fs - 2*Tl;

figure;
subplot(3,2,1);
plot(t, s_rect, 'LineWidth', 2);
xlim([-Tl Tl]); grid on;

subplot(3,2,2);
plot(t, s_tri1, 'LineWidth', 2);
xlim([-Tl Tl]); grid on;

subplot(3,2,3:4);
plot(t_c, c7, 'LineWidth', 2);
xlim([-2*tau 2*T]); grid on;

figure;
subplot(3,2,1);
plot(t, s_rect, 'LineWidth', 2);
xlim([-Tl Tl]); grid on;

subplot(3,2,2);
plot(t, s_tri2, 'LineWidth', 2);
xlim([-Tl Tl]); grid on;

subplot(3,2,3:4);
plot(t_c, c8, 'LineWidth', 2);
xlim([-2*tau 2*T]); grid on;


%2

pkg load signal;

Fs = 200;
A = 1;
tau = 2;
Tl = 10;
t = -Tl:1/Fs:Tl;

s = zeros(size(t));
s(abs(t) <= tau/2) = A;

c = xcorr(s, s) / Fs;

N = 2^13;
freq = (-N/2:N/2-1) * (Fs/N);

S_f = fftshift(fft(s, N)) / Fs;
Es_from_signal = abs(S_f).^2;

C_f = fftshift(fft(c, N)) / Fs;
Es_from_ccf = real(C_f);

figure;
subplot(2,1,1);
plot(freq, Es_from_signal, 'b', 'LineWidth', 2);
xlim([-3 3]);
grid on;

subplot(2,1,2);
plot(freq, Es_from_ccf, 'r--', 'LineWidth', 2);
xlim([-3 3]);
grid on;

max_diff = max(abs(Es_from_signal - Es_from_ccf));
disp(max_diff);
