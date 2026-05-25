clear; clc; close all;

Fs = 1000;
Tau = 2;
t = -Tau:1/Fs:Tau;
A = 1;

% ========== 1. Прямоугольный импульс ==========
tau_rect = 1;
s1 = A * (abs(t) <= tau_rect/2);
[B1, lag] = xcorr(s1, 'normalized');
lag_time = lag/Fs;

% ========== 2. Несимметричный треугольный импульс ==========
T_tri = 1;
s2 = zeros(size(t));
for i = 1:length(t)
    if t(i) >= 0 && t(i) <= T_tri
        s2(i) = A * (t(i)/T_tri);
    end
end
[B2, ~] = xcorr(s2, 'normalized');

% ========== 3. Симметричный треугольный импульс ==========
T_sym = 1;
s3 = A * max(0, 1 - abs(t)/T_sym);
[B3, ~] = xcorr(s3, 'normalized');

% ========== 4. Односторонний экспоненциальный ==========
a = 5;
s4 = A * exp(-a * t) .* (t >= 0);
[B4, ~] = xcorr(s4, 'normalized');

% ========== 5. Двусторонний экспоненциальный ==========
a2 = 5;
s5 = A * exp(-a2 * abs(t));
[B5, ~] = xcorr(s5, 'normalized');

% ========== 6. Меандр (скважность 2) ==========
T_meandr = 1;
t_long = -4:1/Fs:4;
s6_long = A * (mod(t_long + T_meandr/4, T_meandr) < T_meandr/2);
s6 = s6_long(abs(t_long) <= Tau);
[B6, ~] = xcorr(s6, 'normalized');

% ========== ВЫВОД: сигналы ==========
figure('Position', [50, 200, 1400, 900]);
for k = 1:6
    subplot(2,3,k);
    eval(['plot(t, s' num2str(k) ', ''b'', ''LineWidth'', 1.5);']);
    grid on; xlabel('t (сек)'); ylabel('s(t)');
    switch k
        case 1, title('1. Прямоугольный импульс');
        case 2, title('2. Несимметричный треугольный');
        case 3, title('3. Симметричный треугольный');
        case 4, title('4. Односторонний экспоненциальный');
        case 5, title('5. Двусторонний экспоненциальный');
        case 6, title('6. Меандр (скважность 2)');
    end
    xlim([-1.5, 1.5]);
end
sgtitle('Сигналы для корреляционного анализа (ЛР №2, п.3)', 'FontSize', 14);

% ========== ВЫВОД: автокорреляционные функции ==========
figure('Position', [50, 200, 1400, 900]);
for k = 1:6
    subplot(2,3,k);
    eval(['plot(lag_time, B' num2str(k) ', ''r'', ''LineWidth'', 1.5);']);
    grid on; xlabel('\tau (сек)'); ylabel('B(\tau)/max');
    title(['АКФ сигнала ' num2str(k)]);
    xlim([-1.5, 1.5]);
    ylim([-0.2, 1.1]);
end
sgtitle('Автокорреляционные функции (нормированные)', 'FontSize', 14);


%2

clear; clc; close all;

Fs = 2000;
t = -2:1/Fs:2;
A = 1; 
tau = 0.5;

s = A * (abs(t) <= tau/2);

% ========== Прямое вычисление АКФ (БЕЗ нормировки) ==========
[B_direct, lag] = xcorr(s, 'unbiased');
lag_time = lag/Fs;

% ========== Через энергетический спектр ==========
N = length(s);
S = fft(s);
E_spec = abs(S).^2 / N;      % энергетический спектр
B_from_spec = real(ifft(E_spec));
B_from_spec = fftshift(B_from_spec);

% ========== ВИЗУАЛИЗАЦИЯ ==========
figure('Position', [100, 100, 1200, 900]);

subplot(3,2,1);
plot(t, s, 'b', 'LineWidth', 1.5);
grid on; title('1. Прямоугольный импульс s(t)');
xlabel('t (сек)'); ylabel('s(t)'); xlim([-1, 1]);

subplot(3,2,2);
f = (0:N-1)*Fs/N;
plot(f(1:N/2), E_spec(1:N/2), 'g', 'LineWidth', 1.5);
grid on; title('2. Энергетический спектр |S(f)|^2');
xlabel('f (Гц)'); ylabel('|S(f)|^2');
xlim([0, 50]);

subplot(3,2,3);
plot(lag_time, B_direct, 'r', 'LineWidth', 1.5);
grid on; title('3. АКФ (прямое вычисление: xcorr)');
xlabel('\tau (сек)'); ylabel('B(\tau)'); xlim([-1.5, 1.5]);

subplot(3,2,4);
plot(lag_time, B_from_spec, 'm--', 'LineWidth', 1.5);
grid on; title('4. АКФ (через спектр: ifft(|S(f)|^2))');
xlabel('\tau (сек)'); ylabel('B(\tau)'); xlim([-1.5, 1.5]);

subplot(3,2,5);
plot(lag_time, B_direct - B_from_spec', 'k', 'LineWidth', 1);
grid on; title('5. Разница между двумя методами (должна быть ~0)');
xlabel('\tau (сек)'); ylabel('Погрешность'); xlim([-1.5, 1.5]);

subplot(3,2,6);
plot(f(1:N/2), E_spec(1:N/2), 'g', 'LineWidth', 1.5); hold on;
title('6. Энергетический спектр с параметрами');
xlabel('f (Гц)'); ylabel('|S(f)|^2');
xlim([0, 50]);
text(15, max(E_spec(1:N/2))*0.7, sprintf('Ширина главного лепестка: ~%.2f Гц', 1/tau), 'FontSize', 10);
grid on;

sgtitle('Проверка теоремы Винера-Хинчина: АКФ ⟷ Энергетический спектр', 'FontSize', 12);

%3

clear; clc; close all;

Fs = 2000;
t = -2:1/Fs:2;
A = 1; tau = 0.5;

s = A * (abs(t) <= tau/2);

% Прямое вычисление АКФ
[B_direct, lag] = xcorr(s, 'normalized');
lag_time = lag/Fs;

% Через энергетический спектр (теорема Винера-Хинчина)
N = length(s);
S = fft(s);
E_spec = abs(S).^2 / N;      % энергетический спектр
B_from_spec = real(ifft(E_spec));
B_from_spec = fftshift(B_from_spec);
B_from_spec = B_from_spec / max(B_from_spec);

% Обрезаем до той же длины, что и B_direct
B_from_spec = B_from_spec(1:length(B_direct));

figure('Position', [100, 100, 1000, 800]);

subplot(2,2,1);
plot(t, s, 'b', 'LineWidth', 1.5);
grid on; title('Прямоугольный импульс (A=1, \tau=0.5)');
xlabel('t (сек)'); ylabel('s(t)'); xlim([-1, 1]);

subplot(2,2,2);
f = (0:N/2-1)*Fs/N;
plot(f, E_spec(1:N/2), 'g', 'LineWidth', 1.5);
grid on; title('Энергетический спектр |S(f)|^2');
xlabel('f (Гц)'); ylabel('|S(f)|^2');

subplot(2,2,3);
plot(lag_time, B_direct, 'r', 'LineWidth', 1.5);
grid on; title('АКФ (прямое вычисление, xcorr)');
xlabel('\tau (сек)'); ylabel('B(\tau)/max'); xlim([-1.5, 1.5]);

subplot(2,2,4);
plot(lag_time, B_from_spec, 'm--', 'LineWidth', 1.5);
grid on; title('АКФ (через спектр: ifft(|S(f)|^2))');
xlabel('\tau (сек)'); ylabel('B(\tau)/max'); xlim([-1.5, 1.5]);

sgtitle('Проверка связи АКФ и энергетического спектра (теорема Винера-Хинчина)', 'FontSize', 12);