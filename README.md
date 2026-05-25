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
