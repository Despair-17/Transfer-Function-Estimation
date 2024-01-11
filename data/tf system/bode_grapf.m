clear; clc; close all
%%
sysd_axis1 = load('sysd_axis1.mat').sysd_axis;
sysd_axis2 = load('sysd_axis2.mat').sysd_axis;
sysd_axis3 = load('sysd_axis3.mat').sysd_axis;
sysd_axis4 = load('sysd_axis4.mat').sysd_axis;

sysd_rail1 = load('sysd_rail1.mat').sysd_rail;
sysd_rail2 = load('sysd_rail2.mat').sysd_rail;
sysd_rail3 = load('sysd_rail3.mat').sysd_rail;
sysd_rail4 = load('sysd_rail4.mat').sysd_rail;

sysd_diag1 = load('sysd_diag1.mat').sysd_diag;
sysd_diag2 = load('sysd_diag2.mat').sysd_diag;
sysd_diag3 = load('sysd_diag3.mat').sysd_diag;
sysd_diag4 = load('sysd_diag4.mat').sysd_diag;
% 
% figure(10)
% bode(sysd_axis1, sysd_axis2, sysd_axis3, sysd_axis4)
% 
% figure(11)
% bode(sysd_rail1, sysd_rail2, sysd_rail3, sysd_rail4)
% 
% figure(13)
% bode(sysd_diag1, sysd_diag2, sysd_diag3, sysd_diag4)

%% axis
figure(1);

% Задаем передаточные функции
sysd1 = sysd_axis1;
sysd2 = sysd_axis2;
sysd3 = sysd_axis3;
sysd4 = sysd_axis4;

% Получаем амплитуду, фазу и частоту для каждой передаточной функции
[bode_mag1, bode_phase1, freq1] = bode(sysd1);
[bode_mag2, bode_phase2, freq2] = bode(sysd2);
[bode_mag3, bode_phase3, freq3] = bode(sysd3);
[bode_mag4, bode_phase4, freq4] = bode(sysd4);

% Строим ФЧХ
subplot(2,1,2);
semilogx(freq1, squeeze(bode_phase1), freq2, squeeze(bode_phase2), freq3, squeeze(bode_phase3), freq4, squeeze(bode_phase4), 'LineWidth', 1);
grid on;
xlabel('Частота (рад/с)');
ylabel('Фаза (градусы)');
legend('sysd\_axis_1', 'sysd\_axis_2', 'sysd\_axis_3', 'sysd\_axis_4');

% Строим АЧХ
subplot(2,1,1);
semilogx(freq1, 20*log10(squeeze(bode_mag1)), freq2, 20*log10(squeeze(bode_mag2)), freq3, 20*log10(squeeze(bode_mag3)), freq4, 20*log10(squeeze(bode_mag4)), 'LineWidth', 1);
grid on;
xlabel('Частота (рад/с)');
ylabel('Амплитуда (дБ)');
legend('sysd\_axis_1', 'sysd\_axis_2', 'sysd\_axis_3', 'sysd\_axis_4');

%% rail
figure(2);

% Задаем передаточные функции
sysd1 = sysd_rail1;
sysd2 = sysd_rail2;
sysd3 = sysd_rail3;
sysd4 = sysd_rail4;

% Получаем амплитуду, фазу и частоту для каждой передаточной функции
[bode_mag1, bode_phase1, freq1] = bode(sysd1);
[bode_mag2, bode_phase2, freq2] = bode(sysd2);
[bode_mag3, bode_phase3, freq3] = bode(sysd3);
[bode_mag4, bode_phase4, freq4] = bode(sysd4);

% Строим ФЧХ
subplot(2,1,2);
semilogx(freq1, squeeze(bode_phase1), freq2, squeeze(bode_phase2), freq3, squeeze(bode_phase3), freq4, squeeze(bode_phase4), 'LineWidth', 1);
grid on;
xlabel('Частота (рад/с)');
ylabel('Фаза (градусы)');
legend('sysd\_rail_1', 'sysd\_rail_2', 'sysd\_rail_3', 'sysd\_rail_4');

% Строим АЧХ
subplot(2,1,1);
semilogx(freq1, 20*log10(squeeze(bode_mag1)), freq2, 20*log10(squeeze(bode_mag2)), freq3, 20*log10(squeeze(bode_mag3)), freq4, 20*log10(squeeze(bode_mag4)), 'LineWidth', 1);
grid on;
xlabel('Частота (рад/с)');
ylabel('Амплитуда (дБ)');
legend('sysd\_rail_1', 'sysd\_rail_2', 'sysd\_rail_3', 'sysd\_rail_4');

%% diag
figure(3);

% Задаем передаточные функции
sysd1 = sysd_diag1;
sysd2 = sysd_diag2;
sysd3 = sysd_diag3;
sysd4 = sysd_diag4;

% Получаем амплитуду, фазу и частоту для каждой передаточной функции
[bode_mag1, bode_phase1, freq1] = bode(sysd1);
[bode_mag2, bode_phase2, freq2] = bode(sysd2);
[bode_mag3, bode_phase3, freq3] = bode(sysd3);
[bode_mag4, bode_phase4, freq4] = bode(sysd4);

a = 1; b = 1; c = 1; d = 1;
% Строим ФЧХ
subplot(2,1,2);
semilogx(freq1, squeeze(bode_phase1) - a * 360, freq2, squeeze(bode_phase2) - b * 360, freq3, squeeze(bode_phase3) - c * 360, freq4, squeeze(bode_phase4)  - d * 360, 'LineWidth', 1);
grid on;
xlabel('Частота (рад/с)');
ylabel('Фаза (градусы)');
legend('sysd\_diag_1', 'sysd\_diag_2', 'sysd\_diag_3', 'sysd\_diag_4');

% Строим АЧХ
subplot(2,1,1);
semilogx(freq1, 20*log10(squeeze(bode_mag1)), freq2, 20*log10(squeeze(bode_mag2)), freq3, 20*log10(squeeze(bode_mag3)), freq4, 20*log10(squeeze(bode_mag4)), 'LineWidth', 1);
grid on;
xlabel('Частота (рад/с)');
ylabel('Амплитуда (дБ)');
legend('sysd\_diag_1', 'sysd\_diag_2', 'sysd\_diag_3', 'sysd\_diag_4');