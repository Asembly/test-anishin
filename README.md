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

Fs = 100;
dt = 1/Fs;
t = -3:dt:3;

s_rect = zeros(size(t));
s_rect(abs(t) <= 1) = 1;

s_tri = zeros(size(t));
s_tri(abs(t) <= 1) = 1 - abs(t(abs(t) <= 1));

B_rect = xcorr(s_rect, s_rect) * dt;
B_tri = xcorr(s_tri, s_tri) * dt;

t_c = (0:length(B_rect)-1)/Fs - 6;

figure;

subplot(2,2,1);
plot(t, s_rect, 'b', 'LineWidth', 2);
title('1. Прямоугольный импульс s(t)');
xlabel('t, сек'); ylabel('Амплитуда'); grid on;
ylim([-0.2 1.2]); xlim([-2 2]);

subplot(2,2,2);
plot(t_c, B_rect, 'b', 'LineWidth', 2);
title('2. АКФ прямоугольного импульса');
xlabel('\tau, сек'); ylabel('B(\tau)'); grid on;
xlim([-2.5 2.5]);

subplot(2,2,3);
plot(t, s_tri, 'm', 'LineWidth', 2);
title('3. Треугольный импульс s(t)');
xlabel('t, сек'); ylabel('Амплитуда'); grid on;
ylim([-0.2 1.2]); xlim([-2 2]);

subplot(2,2,4);
plot(t_c, B_tri, 'm', 'LineWidth', 2);
title('4. АКФ треугольного импульса');
xlabel('\tau, сек'); ylabel('B(\tau)'); grid on;
xlim([-2.5 2.5]);


// Задание 3
[31.05.2026 20:38] Tylen: % Лабораторная работа №3 (Octave)
% Теорема Котельникова
% Выполните скрипт по частям или целиком.

clear all; close all; clc;

% ========== Задание 1 ==========
% Гармоническое колебание: f1 = 1 Гц, A1 = 1, фаза pi/4, T1 = 10 с
f1 = 1; A1 = 1; phi1 = pi/4; T1 = 10;
% Аналитическая функция (для fplot)
f_analog = @(t) A1 * sin(2*pi*f1*t + phi1);

% Частоты дискретизации
Fs_list = [1, 3, 5, 20];
% Дискретное время для каждого Fs (от 0 до T1, не включая T1, чтобы получить N = T1*Fs отсчётов)
t_cont = linspace(0, T1, 1000); % для отображения аналогового сигнала (плавно)

figure('Name', 'Задание 1: Сигналы');
for i = 1:length(Fs_list)
    Fs = Fs_list(i);
    t = 0 : 1/Fs : T1-1/Fs;   % массив дискретного времени
    s = A1 * sin(2*pi*f1*t + phi1);
    subplot(2,2,i);
    % рисуем аналоговый сигнал красной линией
    plot(t_cont, f_analog(t_cont), 'r-', 'LineWidth', 1.5); hold on;
    % рисуем дискретные отсчёты
    stem(t, s, 'b', 'MarkerSize', 4, 'LineWidth', 1);
    xlabel('t, с'); ylabel('s(t)');
    title(sprintf('Fs = %d Гц', Fs));
    grid on;
    hold off;
end
sgtitle('Дискретизация гармонического сигнала');

% Спектры для каждого Fs
figure('Name', 'Задание 1: Спектры');
for i = 1:length(Fs_list)
    Fs = Fs_list(i);
    t = 0 : 1/Fs : T1-1/Fs;
    s = A1 * sin(2*pi*f1*t + phi1);
    N = length(s);
    Y = fft(s) / N;
    f_axis = (0:N-1) * Fs / N;   % полная ось частот
    % Берём только положительные частоты до Fs/2
    half = floor(N/2) + 1;
    Y_half = Y(1:half);
    f_half = f_axis(1:half);
    % Амплитудный спектр (умножаем на 2, кроме постоянной составляющей)
    amp = 2 * abs(Y_half);
    amp(1) = amp(1) / 2;   % поправка для нулевой частоты (если есть)
    subplot(2,2,i);
    stem(f_half, amp, 'filled', 'MarkerSize', 4);
    xlabel('f, Гц'); ylabel('|S(f)|');
    title(sprintf('Спектр, Fs = %d Гц', Fs));
    xlim([0, Fs/2 + 0.5]);
    grid on;
end
sgtitle('Амплитудные спектры дискретных сигналов');

% ========== Задание 2 ==========
% Сигнал: сумма двух гармоник: f1=1 Гц, f2=3 Гц, A1=0.7, A2=0.9, фазы 0 и pi/8
f1_2 = 1; f2_2 = 3; A1_2 = 0.7; A2_2 = 0.9; phi1_2 = 0; phi2_2 = pi/8;
f_sum = @(t) A1_2 * sin(2*pi*f1_2*t + phi1_2) + A2_2 * sin(2*pi*f2_2*t + phi2_2);

% Те же частоты дискретизации
Fs_list2 = [1, 3, 5, 20];
T1_2 = 10;   % длительность
t_cont2 = linspace(0, T1_2, 1000);

figure('Name', 'Задание 2: Сигналы (сумма гармоник)');
for i = 1:length(Fs_list2)
    Fs = Fs_list2(i);
    t = 0 : 1/Fs : T1_2-1/Fs;
    s = f_sum(t);
    subplot(2,2,i);
    plot(t_cont2, f_sum(t_cont2), 'r-', 'LineWidth', 1.5); hold on;
    stem(t, s, 'b', 'MarkerSize', 4);
    xlabel('t, с'); ylabel('s(t)');
    title(sprintf('Fs = %d Гц', Fs));
    grid on;
    hold off;
end
sgtitle('Дискретизация суммы двух гармоник');

% Спектры для задания 2
figure('Name', 'Задание 2: Спектры');
for i = 1:length(Fs_list2)
    Fs = Fs_list2(i);
    t = 0 : 1/Fs : T1_2-1/Fs;
    s = f_sum(t);
    N = length(s);
    Y = fft(s) / N;
    f_axis = (0:N-1) * Fs / N;
    half = floor(N/2) + 1;
    Y_half = Y(1:half);
    f_half = f_axis(1:half);
    amp = 2 * abs(Y_half);
    amp(1) = amp(1) / 2;
    subplot(2,2,i);
    stem(f_half, amp, 'filled', 'MarkerSize', 4);
    xlabel('f, Гц'); ylabel('|S(f)|');
    title(sprintf('Спектр, Fs = %d Гц', Fs));
    xlim([0, Fs/2 + 1]);
    grid on;
end
sgtitle('Амплитудные спектры (сумма гармоник)');

% ========== Задание 3 ==========
% Две гармоники с девиацией частоты Δf, начальные фазы 0 и pi/8
% Частоты: f0 = 5 Гц (некоторая центральная), и f0+Δf
% По заданию: "сумма 2-х гармонических колебаний с девиацией частоты Δf"
% По умолчанию возьмём центральную частоту 5 Гц, тогда частоты f_c - Δf/2 и f_c + Δf/2
% Но в тексте просто указана Δf – будем считать, что частоты f1 = f0, f2 = f0+Δf.
% Для определённости f0 = 5 Гц.

f0 = 5;  % центральная частота (базовая)
A = 1;   % амплитуды одинаковы
phi = [0, pi/8];
Delta_f_list = [0.5, 0.3, 0.1, 0.01];
[31.05.2026 20:38] Tylen: % Для каждой Δf подбираем T1 так, чтобы разрешение БПФ позволяло различить пики.
% Разрешение по частоте: df = 1 / T1. Чтобы различать частоты, нужно T1 >= 1/Δf.
% Возьмём T1 = 2/Δf для надёжности, но с минимальным значением 10 с.
fprintf('\nЗадание 3: Определение длительности сигнала для разрешения спектральных составляющих\n');
for k = 1:length(Delta_f_list)
    df = Delta_f_list(k);
    T1_rec = max(10, 2/df);   % рекомендуемая длительность
    fprintf('Δf = %.3f Гц -> рекомендуемая длительность сигнала T1 >= %.2f с\n', df, T1_rec);
end

% Для демонстрации выберем T1 = 2/Δf (но не менее 10) и построим спектры.
figure('Name', 'Задание 3: Спектры для разных Δf');
for k = 1:length(Delta_f_list)
    df = Delta_f_list(k);
    T1_cur = max(10, 2/df);
    % Частоты: f0 и f0+df (можно и f0-df/2, f0+df/2 – суть та же)
    fA = f0;
    fB = f0 + df;
    % Аналитическая функция
    f_sig = @(t) A * sin(2*pi*fA*t + phi(1)) + A * sin(2*pi*fB*t + phi(2));
    % Частота дискретизации: возьмём заведомо больше 2*(f0+df)
    Fs = 100;   % высокая, чтобы не было наложения
    t = 0 : 1/Fs : T1_cur - 1/Fs;
    s = f_sig(t);
    N = length(s);
    Y = fft(s) / N;
    f_axis = (0:N-1) * Fs / N;
    half = floor(N/2) + 1;
    Y_half = Y(1:half);
    f_half = f_axis(1:half);
    amp = 2 * abs(Y_half);
    amp(1) = amp(1) / 2;
    subplot(2,2,k);
    plot(f_half, amp, 'b-', 'LineWidth', 1.5);
    xlabel('f, Гц'); ylabel('|S(f)|');
    title(sprintf('Δf = %.3f Гц, T1 = %.1f с', df, T1_cur));
    xlim([f0-2, f0+df+2]);
    grid on;
    % отметим теоретические частоты
    hold on;
    plot([fA fA], ylim, 'r--', 'LineWidth', 1);
    plot([fB fB], ylim, 'r--', 'LineWidth', 1);
    hold off;
end
sgtitle('Разрешение двух близких частот при разной девиации');

//Задание 4
[31.05.2026 20:40] Tylen: % Лабораторная работа №4 (Octave)
% Амплитудная, частотная и фазовая модуляция
clear all; close all; clc;

% Общие параметры
Fs = 1000;          % частота дискретизации (Гц)
T_total = 2;        % общая длительность сигнала (с)
t = 0 : 1/Fs : T_total - 1/Fs;   % временная сетка
% Несущее колебание
fc = 50;            % несущая частота (Гц)
Ac = 1;             % амплитуда несущей
carrier = Ac * cos(2*pi*fc*t);

% Коэффициенты модуляции (подбираются для наглядности)
k_AM = 0.8;         % глубина АМ (чтобы огибающая не переходила через ноль)
k_FM = 50;          % девиация частоты для ЧМ (Гц/ед. амплитуды модулирующего)
k_PM = pi/2;        % девиация фазы для ФМ (рад/ед. амплитуды)

% ========== 1. Несимметричный треугольный импульс ==========
% Параметры: A = 1, T = 1 с (длительность импульса)
A_tri = 1;
T_tri = 1;
% Формируем модулирующий сигнал s_mod (треугольник на интервале [0, T_tri])
s_mod_tri = zeros(size(t));
idx = (t >= 0) & (t <= T_tri);
s_mod_tri(idx) = A_tri * (t(idx) / T_tri);   % возрастающая прямая

% АМ для треугольного импульса
s_AM_tri = (1 + k_AM * s_mod_tri) .* carrier;
% ЧМ: интеграл от модулирующего сигнала
integral_tri = cumtrapz(t, s_mod_tri);   % численное интегрирование
phi_FM_tri = 2*pi*k_FM * integral_tri;
s_FM_tri = Ac * cos(2*pi*fc*t + phi_FM_tri);
% ФМ: фаза пропорциональна модулирующему сигналу
phi_PM_tri = k_PM * s_mod_tri;
s_PM_tri = Ac * cos(2*pi*fc*t + phi_PM_tri);

% Визуализация для треугольного импульса
figure('Name', 'Треугольный импульс: модуляция');
% Осциллограммы (первые 0.5 с для наглядности)
t_zoom = t(t <= 0.5);
s_mod_zoom = s_mod_tri(t <= 0.5);
carrier_zoom = carrier(t <= 0.5);
s_AM_zoom = s_AM_tri(t <= 0.5);
s_FM_zoom = s_FM_tri(t <= 0.5);
s_PM_zoom = s_PM_tri(t <= 0.5);

subplot(3,3,1); plot(t_zoom, s_mod_zoom, 'b'); title('Модулирующий сигнал');
xlabel('t, с'); ylabel('s_m'); grid on;
subplot(3,3,2); plot(t_zoom, carrier_zoom, 'k'); title('Несущая');
xlabel('t, с'); ylabel('Несущая'); grid on;
subplot(3,3,3); plot(t_zoom, s_AM_zoom, 'r'); title('АМ-сигнал');
xlabel('t, с'); ylabel('s_{AM}'); grid on;
subplot(3,3,4); plot(t_zoom, s_FM_zoom, 'g'); title('ЧМ-сигнал');
xlabel('t, с'); ylabel('s_{FM}'); grid on;
subplot(3,3,5); plot(t_zoom, s_PM_zoom, 'm'); title('ФМ-сигнал');
xlabel('t, с'); ylabel('s_{PM}'); grid on;
% Спектры (полные)
compute_and_plot_spectrum(s_mod_tri, Fs, 'Спектр модулирующего', 3,6);
compute_and_plot_spectrum(s_AM_tri, Fs, 'Спектр АМ', 3,7);
compute_and_plot_spectrum(s_FM_tri, Fs, 'Спектр ЧМ', 3,8);
compute_and_plot_spectrum(s_PM_tri, Fs, 'Спектр ФМ', 3,9);
sgtitle('Модуляция треугольным импульсом');

% ========== 2. Одиночный прямоугольный импульс (длительности τ = 0.1, 0.5, 1 с) ==========
tau_list = [0.1, 0.5, 1];
for idx_tau = 1:length(tau_list)
    tau = tau_list(idx_tau);
    s_mod_rect = zeros(size(t));
    idx_rect = (t >= 0) & (t <= tau);
    s_mod_rect(idx_rect) = 1;   % амплитуда 1
    
    % АМ
    s_AM_rect = (1 + k_AM * s_mod_rect) .* carrier;
    % ЧМ
    integral_rect = cumtrapz(t, s_mod_rect);
    phi_FM_rect = 2*pi*k_FM * integral_rect;
    s_FM_rect = Ac * cos(2*pi*fc*t + phi_FM_rect);
    % ФМ
    phi_PM_rect = k_PM * s_mod_rect;
    s_PM_rect = Ac * cos(2*pi*fc*t + phi_PM_rect);
    
    figure('Name', sprintf('Прямоугольный импульс, τ = %.1f с', tau));
    % Осциллограммы (первые 0.5 с или до 2*tau)
    t_zoom_rect = t(t <= min(0.8, 2*tau));
    s_mod_zoom_rect = s_mod_rect(t <= min(0.8, 2*tau));
    carrier_zoom_rect = carrier(t <= min(0.8, 2*tau));
    s_AM_zoom_rect = s_AM_rect(t <= min(0.8, 2*tau));
    s_FM_zoom_rect = s_FM_rect(t <= min(0.8, 2*tau));
    s_PM_zoom_rect = s_PM_rect(t <= min(0.8, 2*tau));
    
    subplot(3,3,1); plot(t_zoom_rect, s_mod_zoom_rect, 'b'); title('Модулирующий');
    xlabel('t, с'); grid on;
    subplot(3,3,2); plot(t_zoom_rect, carrier_zoom_rect, 'k'); title('Несущая');
[31.05.2026 20:40] Tylen: xlabel('t, с'); grid on;
    subplot(3,3,3); plot(t_zoom_rect, s_AM_zoom_rect, 'r'); title('АМ');
    xlabel('t, с'); grid on;
    subplot(3,3,4); plot(t_zoom_rect, s_FM_zoom_rect, 'g'); title('ЧМ');
    xlabel('t, с'); grid on;
    subplot(3,3,5); plot(t_zoom_rect, s_PM_zoom_rect, 'm'); title('ФМ');
    xlabel('t, с'); grid on;
    % Спектры
    compute_and_plot_spectrum(s_mod_rect, Fs, 'Спектр модулирующего', 3,6);
    compute_and_plot_spectrum(s_AM_rect, Fs, 'Спектр АМ', 3,7);
    compute_and_plot_spectrum(s_FM_rect, Fs, 'Спектр ЧМ', 3,8);
    compute_and_plot_spectrum(s_PM_rect, Fs, 'Спектр ФМ', 3,9);
    sgtitle(sprintf('Прямоугольный импульс, длительность τ = %.1f с', tau));
end

% ========== 3. Гармонический модулирующий сигнал (частота Fm = 1 Гц) ==========
Fm = 1;
s_mod_sin = sin(2*pi*Fm*t);   % амплитуда 1
% Для АМ нужна однополярная огибающая, поэтому добавим постоянную составляющую 1
% s_AM_sin = (1 + k_AM * s_mod_sin) .* carrier
% Но для синусоидальной модуляции это стандартная АМ.
s_AM_sin = (1 + k_AM * s_mod_sin) .* carrier;
% ЧМ: интеграл от синуса даёт -косинус, но используем численное интегрирование
integral_sin = cumtrapz(t, s_mod_sin);
phi_FM_sin = 2*pi*k_FM * integral_sin;
s_FM_sin = Ac * cos(2*pi*fc*t + phi_FM_sin);
% ФМ
phi_PM_sin = k_PM * s_mod_sin;
s_PM_sin = Ac * cos(2*pi*fc*t + phi_PM_sin);

figure('Name', 'Гармонический модулирующий сигнал (Fm=1 Гц)');
t_zoom_sin = t(t <= 0.5);
s_mod_sin_zoom = s_mod_sin(t <= 0.5);
carrier_zoom_sin = carrier(t <= 0.5);
s_AM_sin_zoom = s_AM_sin(t <= 0.5);
s_FM_sin_zoom = s_FM_sin(t <= 0.5);
s_PM_sin_zoom = s_PM_sin(t <= 0.5);

subplot(3,3,1); plot(t_zoom_sin, s_mod_sin_zoom, 'b'); title('Модулирующий');
xlabel('t, с'); grid on;
subplot(3,3,2); plot(t_zoom_sin, carrier_zoom_sin, 'k'); title('Несущая');
xlabel('t, с'); grid on;
subplot(3,3,3); plot(t_zoom_sin, s_AM_sin_zoom, 'r'); title('АМ');
xlabel('t, с'); grid on;
subplot(3,3,4); plot(t_zoom_sin, s_FM_sin_zoom, 'g'); title('ЧМ');
xlabel('t, с'); grid on;
subplot(3,3,5); plot(t_zoom_sin, s_PM_sin_zoom, 'm'); title('ФМ');
xlabel('t, с'); grid on;
% Спектры
compute_and_plot_spectrum(s_mod_sin, Fs, 'Спектр модулирующего', 3,6);
compute_and_plot_spectrum(s_AM_sin, Fs, 'Спектр АМ', 3,7);
compute_and_plot_spectrum(s_FM_sin, Fs, 'Спектр ЧМ', 3,8);
compute_and_plot_spectrum(s_PM_sin, Fs, 'Спектр ФМ', 3,9);
sgtitle('Модуляция гармоническим сигналом 1 Гц');

% Вспомогательная функция для расчёта и построения амплитудного спектра
function compute_and_plot_spectrum(signal, Fs, title_str, sub_row, sub_col)
    N = length(signal);
    Y = fft(signal) / N;
    f_axis = (0:N-1) * Fs / N;
    half = floor(N/2) + 1;
    Y_half = Y(1:half);
    f_half = f_axis(1:half);
    amp = 2 * abs(Y_half);
    amp(1) = amp(1) / 2;   % поправка постоянной составляющей
    subplot(sub_row, sub_col, sub_col * (sub_row-1) + sub_col); % простая индексация
    % Используем plot вместо stem для лучшего восприятия при большом числе частот
    plot(f_half, amp, 'b-', 'LineWidth', 1);
    xlabel('f, Гц'); ylabel('|S(f)|');
    title(title_str);
    xlim([0, Fs/2]);
    grid on;
end

Индивидуальное задание 2
x =[1, 0, 1, 2, 1, 0, 1];
y =[1, 1, 1];

% 2. Вычисляем взаимную корреляцию
R = conv(x, fliplr(y));

% 3. Создаем ось сдвигов (m)
m = -(length(y)-1) : (length(x)-1);

% 4. Вывод значений в консоль
fprintf('Значения ВКФ:\n');
fprintf('%d  ', R);
fprintf('\n');

% 5. Построение графика
figure(1);
stem(m, R, 'filled', 'LineWidth', 2, 'MarkerSize', 8);
grid on;
xlabel('Сдвиг (m)');
ylabel('R_{xy}[m]');
title('Взаимная корреляционная функция сигналов');
xlim([min(m)-1, max(m)+1]);
ylim([0, max(R)+1]);

Индивидуальное задание 3
