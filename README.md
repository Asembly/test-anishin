clear; clc; close all;

T = 1.0;
t = linspace(-0.2, 2.8, 4000);
s = zeros(size(t));

for i = 1:length(t)
    t_mod = mod(t(i), T);
    
    if t_mod >= 0 && t_mod < (T / 3)
        w_neg = T / 3;
        s(i) = -1.0 + abs(t_mod - w_neg/2) / (w_neg/2);
    else
        if t_mod < (3/6)
            s(i) = 1.2 * (t_mod - 2/6) / (1/6);
        elseif t_mod < (4/6)
            s(i) = 1.2 - 0.6 * (t_mod - 3/6) / (1/6);
        elseif t_mod < (5/6)
            s(i) = 0.6 + 0.6 * (t_mod - 4/6) / (1/6);
        else
            s(i) = 1.2 * (1.0 - t_mod) / (1/6);
        end
    end
end

Fs = 1 / (t(2) - t(1));
N = length(t);
f = (0:N-1) * (Fs / N);
S_fft = abs(fft(s)) / N;

max_freq_idx = find(f > 25, 1);

figure;

subplot(2,1,1);
plot(t, s, 'b', 'LineWidth', 1.2);
hold on;
plot([min(t) max(t)], [0 0], 'k', 'LineWidth', 0.8);
hold off;
grid on;
title('Временной сигнал s(t) (График №14 без зазоров)');
xlabel('Время t');
ylabel('s(t)');
axis([min(t) max(t) -1.3 1.5]);

subplot(2,1,2);
plot(f(1:max_freq_idx), S_fft(1:max_freq_idx), 'r', 'LineWidth', 1.5);
grid on;
title('Амплитудный спектр непрерывного сигнала S(f)');
xlabel('Частота f');
ylabel('|S(f)|');