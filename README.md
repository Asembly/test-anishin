clear; clc; close all;

Fs = 1000;
Tau = 2;
t = -Tau:1/Fs:Tau;

% 1. Прямоугольный импульс
A = 1; tau = 1;
s1 = A * (abs(t) <= tau/2);
B1 = xcorr(s1, 'normalized');
lag = (-length(s1)+1:length(s1)-1)/Fs;

figure('Position', [100, 100, 1200, 800]);
for k = 1:6
    subplot(3,2,k);
    switch k
        case 1  % Прямоугольный
            plot(t, s1, 'b', 'LineWidth', 1.5);
            title('Прямоугольный импульс'); grid on;
        case 2  % Несимметричный треугольный
            T_tri = 1;
            s2 = max(0, A * (1 - t/T_tri)) .* (t >= 0);
            B2 = xcorr(s2, 'normalized');
            plot(t, s2, 'b', 'LineWidth', 1.5);
            title('Несимметричный треугольный'); grid on;
        case 3  % Симметричный треугольный
            T_tri_sym = 1;
            s3 = A * max(0, 1 - abs(t)/T_tri_sym);
            B3 = xcorr(s3, 'normalized');
            plot(t, s3, 'b', 'LineWidth', 1.5);
            title('Симметричный треугольный'); grid on;
        case 4  % Односторонний экспоненциальный
            a = 5;
            s4 = (t >= 0) .* A * exp(-a*t);
            B4 = xcorr(s4, 'normalized');
            plot(t, s4, 'b', 'LineWidth', 1.5);
            title('Односторонний экспоненциальный'); grid on;
        case 5  % Двусторонний экспоненциальный
            a = 5;
            s5 = A * exp(-a*abs(t));
            B5 = xcorr(s5, 'normalized');
            plot(t, s5, 'b', 'LineWidth', 1.5);
            title('Двусторонний экспоненциальный'); grid on;
        case 6  % Меандр
            A_meandr = 1; T_meandr = 1; N_periods = 4;
            t_meandr = -Tau:1/Fs:Tau;
            s6 = A_meandr * square(2*pi/T_meandr * t_meandr, 50);
            B6 = xcorr(s6, 'normalized');
            plot(t_meandr, s6, 'b', 'LineWidth', 1.5);
            title('Меандр (скважность 2)'); grid on;
    end
    xlabel('t (сек)'); ylabel('s(t)');
    xlim([-2, 2]);
end
sgtitle('Сигналы для корреляционного анализа');

% Графики АКФ
figure('Position', [100, 100, 1200, 800]);
for k = 1:6
    subplot(3,2,k);
    switch k
        case 1, B = B1; case 2, B = B2; case 3, B = B3;
        case 4, B = B4; case 5, B = B5; case 6, B = B6;
    end
    plot(lag, B, 'r', 'LineWidth', 1.5);
    title(['АКФ ', num2str(k)]); grid on;
    xlabel('\tau (сек)'); ylabel('B(\tau)/max');
    xlim([-2, 2]);
end
sgtitle('Автокорреляционные функции');


%2
clear; clc; close all;

Fs = 1000;
t = -2:1/Fs:2;

% Сигналы из пункта 7 (прямоугольный и треугольный)
s_rect = (abs(t) <= 0.5);
s_tri = max(0, 1 - abs(t)/1);

% Взаимная корреляция
[B12, lag] = xcorr(s_rect, s_tri, 'normalized');
lag_time = lag/Fs;

figure;
subplot(3,1,1);
plot(t, s_rect, 'b', 'LineWidth', 1.5); hold on;
plot(t, s_tri, 'r', 'LineWidth', 1.5);
legend('Прямоугольный', 'Треугольный');
title('Сигналы s1(t) и s2(t)'); grid on;

subplot(3,1,2);
plot(lag_time, B12, 'g', 'LineWidth', 1.5);
title('Взаимная корреляционная функция B_{12}(\tau)');
xlabel('\tau (сек)'); ylabel('B_{12}(\tau)'); grid on;

subplot(3,1,3);
[B21, lag2] = xcorr(s_tri, s_rect, 'normalized');
plot(lag2/Fs, B21, 'm', 'LineWidth', 1.5);
title('Взаимная корреляционная функция B_{21}(\tau)');
xlabel('\tau (сек)'); ylabel('B_{21}(\tau)'); grid on;

%3
clear; clc; close all;

Fs = 2000;
t = -2:1/Fs:2;
A = 1; tau = 0.5;

s = A * (abs(t) <= tau/2);

% АКФ через xcorr
B_direct = xcorr(s, 'normalized');
lag = (-length(s)+1:length(s)-1)/Fs;

% Энергетический спектр и обратное ПФ
N = length(s);
S = fft(s);
E_spec = abs(S).^2 / N;
B_from_spec = real(ifft(E_spec));
B_from_spec = fftshift(B_from_spec);
B_from_spec = B_from_spec / max(B_from_spec);

figure;
subplot(2,2,1);
plot(t, s, 'b', 'LineWidth', 1.5);
title('Прямоугольный импульс'); grid on;

subplot(2,2,2);
f = (0:N-1)*Fs/N;
plot(f(1:N/2), E_spec(1:N/2), 'r', 'LineWidth', 1.5);
title('Энергетический спектр |S(f)|^2'); grid on;
xlabel('f (Гц)');

subplot(2,2,3);
plot(lag, B_direct, 'g', 'LineWidth', 1.5);
title('АКФ (прямое вычисление)'); grid on;
xlabel('\tau (сек)');

subplot(2,2,4);
plot(lag, B_from_spec, 'm--', 'LineWidth', 1.5);
title('АКФ (через спектр → обратное ПФ)'); grid on;
xlabel('\tau (сек)');