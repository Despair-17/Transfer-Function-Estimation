clear; close all; clc
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Путь к файлам
folder_path = pwd();
main_path = fileparts(mfilename('fullpath'));

% Загурзка данных ручных измерений
Manual_measure = load(fullfile(main_path, "manual_measurement.txt"));
[name_area, ~, ~] = fileparts(folder_path);
[~, idx_area, ~] = fileparts(name_area);

% Загрузка исходных данных и скорректированных данных
[accel, V1mean, V2mean, Ts]= loadAndPreprocess(folder_path, main_path, 0);
load(fullfile(folder_path, "accel_corrected.mat"));

% Нормирование скорректированных данных
accel(:, 1:2) = accel(:,1:2) * (V1mean^2);
accel(:, 3:4) = accel(:,3:4) * (V2mean^2);

accel_corrected(:, 1:2) = accel_corrected(:,1:2);
accel_corrected(:, 3:4) = accel_corrected(:,3:4);

n = 8;
accel(:, 1) = smoothdata(accel(:, 1), 'movmean', n);
accel(:, 2) = smoothdata(accel(:, 2), 'movmean', n);
accel(:, 3) = smoothdata(accel(:, 3), 'movmean', n);
accel(:, 4) = smoothdata(accel(:, 4), 'movmean', n);

% Огибающие сигнала для удобного вычисления максимума
[~, idx] = max(envelope(abs([accel(:, 1), accel(:, 2), ...
                        accel(:, 3), accel(:, 4)]), 100, 'peak'), [], 1);

for i = 1:numel(idx)
    if idx(i) < 501
        idx(i) = idx(i) + 1000;
    end 

    if idx(i) > (length(accel) - 1500)
        idx(i) = idx(i) - 1500;
    end
end

% Диапазон исследуемых участков
len_before = 500; 
len_after = 1499;
range_signal1 = idx(1) - len_before: idx(1) + len_after;
range_signal2 = idx(2) - len_before: idx(2) + len_after;
range_signal3 = idx(3) - len_before: idx(3) + len_after;
range_signal4 = idx(4) - len_before: idx(4) + len_after;

% Выделение участков с полезным сигналом
accel_comb = [accel(range_signal1, 1), ...
              accel(range_signal2, 2), ...
              accel(range_signal3, 3), ...
              accel(range_signal4, 4)];

[idx_corr2, idx_corr4] = index_correction(accel_comb(:,1), ...
                                          accel_comb(:,2), ...
                                          accel_comb(:,3), ...
                                          accel_comb(:,4));
if abs(idx_corr2) > 300
    idx_corr2 = 0;
end

if abs(idx_corr4) > 300
    idx_corr4 = 0;
end

% Совмещение данных проезда одного дефекта впереди и следом идущем колесом
range_signal2 = idx(2) - len_before + idx_corr2: ...
                idx(2) + len_after + idx_corr2;
range_signal4 = idx(4) - len_before + idx_corr4: ...
                idx(4) + len_after + idx_corr4;

accel_comb = [accel(range_signal1, 1), ...
              accel(range_signal2, 2), ...
              accel(range_signal3, 3), ...
              accel(range_signal4, 4)];

accel_corrected_comb = [accel_corrected(range_signal1, 1), ...
                        accel_corrected(range_signal2, 2), ...
                        accel_corrected(range_signal3, 3), ...
                        accel_corrected(range_signal4, 4)];
m = 400;
n = 100;

% Глубина дефекта исходная в мм
h1 = (10^3) * cumsum(cumsum(accel(idx(1) - m:idx(1) + n, 1) * Ts) * Ts);
h2 = (10^3) * cumsum(cumsum(accel(idx(2) - m:idx(2) + n, 2) * Ts) * Ts);
h3 = (10^3) * cumsum(cumsum(accel(idx(3) - m:idx(3) + n, 3) * Ts) * Ts);
h4 = (10^3) * cumsum(cumsum(accel(idx(4) - m:idx(4) + n, 4) * Ts) * Ts);

% Глубина дефекта скорректированная в мм
hcorr1 = (10^3) *cumsum(cumsum(accel_corrected(idx(1) - m:idx(1) + n, 1)...
                                                      * Ts) * Ts);
hcorr2 = (10^3) *cumsum(cumsum(accel_corrected(idx(2) - m:idx(2) + n, 2)...
                                                      * Ts) * Ts);
hcorr3 = (10^3) *cumsum(cumsum(accel_corrected(idx(3) - m:idx(3) + n, 3)...
                                                       * Ts) * Ts);
hcorr4 = (10^3) *cumsum(cumsum(accel_corrected(idx(4) - m:idx(4) + n, 4)...
                                                       * Ts) * Ts);
idx_area = str2double(idx_area) + 1;

image_data  = imread(fullfile(main_path, 'wheels.PNG'));

figure
imshow(image_data);

figure
plot(1:length(accel(:,1)), accel(:,1),'k', 1:length(accel(:,2)), accel(:,2), 'r', ...
    1:length(accel(:,3)), accel(:,3), 'b', 1:length(accel(:,4)), accel(:,4),'g')
    xlabel('Отсчеты')
%     ylabel('Нормированные ускорения, м^-1')
    ylabel('Ускорения, м/с^2')
    legend('1 Следом идущее','2 Впереди идущее', ...
        '3 Следом идущее','4 Впереди идущее')
        legend('1 Впереди идущее','2  Следом идущее', ...
        '3 Впереди идущее','4 Следом идущее')
    grid on
    xlim([0 length(accel(:,1))])
    ylim([min(min(accel(:,1)) - 80, min(accel(:,3)) - 80) ...
        max(max(accel(:,1)) + 80, max(accel(:,3)) + 80)])

figure
plot(h1, 'k')
hold on 
plot(hcorr1, 'b')
xlims = xlim;
line([xlims(1) xlims(2)], [-Manual_measure(idx_area, 1) -Manual_measure(idx_area, 1)], 'Color', 'r'); 
hold off 
grid on
ylabel('Глубина дефекта, мм')
xlabel('Отсчеты')
legend('Исходные данные', 'Скорректированные данные', 'Ручное измерение')
% title('Колесо 1')
ylim([-2 1])

figure
plot(h2, 'k')
hold on 
plot(hcorr2, 'b')
xlims = xlim;
line([xlims(1) xlims(2)], [-Manual_measure(idx_area, 1) -Manual_measure(idx_area, 1)], 'Color', 'r'); 
hold off
grid on
ylabel('Глубина дефекта, мм')
xlabel('Отсчеты')
legend('Исходные данные', 'Скорректированные данные', 'Ручное измерение')
% title('Колесо 2')
ylim([-2 1])

figure
plot(h3, 'k')
hold on 
plot(hcorr3, 'b')
xlims = xlim;
line([xlims(1) xlims(2)], [-Manual_measure(idx_area, 1) -Manual_measure(idx_area, 1)], 'Color', 'r'); 
hold off
grid on
ylabel('Глубина дефекта, мм')
xlabel('Отсчеты')
legend('Исходные данные', 'Скорректированные данные', 'Ручное измерение')
% title('Колесо 3')
ylim([-2 1])

figure
plot(h4, 'k')
hold on 
plot(hcorr4, 'b')
xlims = xlim;
line([xlims(1) xlims(2)], [-Manual_measure(idx_area, 1) -Manual_measure(idx_area, 1)], 'Color', 'r'); 
hold off
grid on
ylabel('Глубина дефекта, мм')
xlabel('Отсчеты')
legend('Исходные данные', 'Скорректированные данные', 'Ручное измерение')
% title('Колесо 4')
ylim([-2 1])

% clear range_signal1 range_signal2 range_signal3 range_signal4 ...
%     idx_corr2 idx_corr4 xlims m n
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Графики 

% Сдвиг сигналов
shift = 300;

% Корреляция исходных
[R12,f1] = xcorr(accel_comb(:, 1), accel_comb(:, 2), shift, 'normalized');
[R34,f2] = xcorr(accel_comb(:, 3), accel_comb(:, 4), shift, 'normalized');

% Корреляция скорректированных
[R12_corr, f1_corr] = xcorr(accel_corrected_comb(:, 1), ...
                      accel_corrected_comb(:, 2), shift, 'normalized');
[R34_corr, f2_corr] = xcorr(accel_corrected_comb(:, 3), ...
                      accel_corrected_comb(:, 4), shift, 'normalized');

% Графики исходных
figure
plot(1:length(accel_comb(:, 1)), accel_comb(:, 1),'k', ...
    1:length(accel_comb(:, 2)), accel_comb(:, 2), 'r')
    xlabel('Отсчеты')
    ylabel('Нормированные ускорения, м^-1')
%     ylabel('Ускорения, м/с^2')
    legend('1 Следом идущее L','2 Впереди идущее L')
    grid on  
    ylim([min(min(accel_comb(:, 1)) - 1, min(accel_comb(:, 3)) - 1) ...
        max(max(accel_comb(:, 1)) + 1, max(accel_comb(:, 3)) + 1)])

figure
plot(1:length(accel_comb(:, 3)), accel_comb(:, 3), 'b', ...
    1:length(accel_comb(:, 4)), accel_comb(:, 4),'g')
    xlabel('Отсчеты')
    ylabel('Нормированные ускорения, м^-1')
%     ylabel('Ускорения, м/с^2')
    legend('3 Следом идущее R','4 Впереди идущее R')
    grid on  
    ylim([min(min(accel_comb(:, 1)) - 1, min(accel_comb(:, 3)) - 1) ...
        max(max(accel_comb(:, 1)) + 1, max(accel_comb(:, 3)) + 1)])

figure
subplot(2, 1, 1)
plot(f1, R12)
    title("До учета")
    xlabel('Отсчеты')
    ylabel('Коэффициент корреляции')
    legend('1 и 2 колесо')
    ylim([-0.7, 1])
    grid on
    set(gca, 'YTick',-0.5:0.2:1)
    

subplot(2, 1, 2)
plot(f2, R34)
    xlabel('Отсчеты')
    ylabel('Коэффициент корреляции')
    legend('3 и 4 колесо')
    ylim([-0.7, 1])
    grid on
    set(gca, 'YTick',-0.5:0.2:1)

% Графики скорректированных

figure
plot(1:length(accel_corrected_comb(:, 1)), accel_corrected_comb(:, 1),'k', ...
     1:length(accel_corrected_comb(:, 2)), accel_corrected_comb(:, 2), 'r')
    xlabel('Отсчеты')
    ylabel('Нормированные ускорения, м^-1')
%     ylabel('Ускорения, м/с^2')
    legend('1 Следом идущее L','2 Впереди идущее L')
    grid on  
    ylim([min(min(accel_comb(:, 1)) - 1, min(accel_comb(:, 3)) - 1) ...
        max(max(accel_comb(:, 1)) + 1, max(accel_comb(:, 3)) + 1)])

figure
plot(1:length(accel_corrected_comb(:, 3)), accel_corrected_comb(:, 3), 'b', ...
     1:length(accel_corrected_comb(:, 4)), accel_corrected_comb(:, 4),'g')
    xlabel('Отсчеты')
    ylabel('Нормированные ускорения, м^-1')
%     ylabel('Ускорения, м/с^2')
    legend('3 Следом идущее R','4 Впереди идущее R')
    grid on  
    ylim([min(min(accel_comb(:, 1)) - 1, min(accel_comb(:, 3)) - 1) ...
        max(max(accel_comb(:, 1)) + 1, max(accel_comb(:, 3)) + 1)])

figure
subplot(2, 1, 1)
plot(f1_corr, R12_corr)
    title("До учета")
    xlabel('Отсчеты')
    ylabel('Коэффициент корреляции')
    legend('1 и 2 колесо')
    ylim([-0.7, 1])
    grid on
    set(gca, 'YTick',-0.5:0.2:1)
    
subplot(2, 1, 2)
plot(f2_corr, R34_corr)
    xlabel('Отсчеты')
    ylabel('Коэффициент корреляции')
    legend('3 и 4 колесо')
    ylim([-0.7, 1])
    grid on
    set(gca, 'YTick',-0.5:0.2:1)

disp("Корреляция до учета...")
disp(['Максимум корреляции сигналов с 1 и 2 датчика: ' ...
    num2str(max(R12))]);
disp(['Максимум корреляции сигналов с 3 и 4 датчика: ' ...
    num2str(max(R34))]);
disp("Корреляция после учета...")
disp(['Максимум корреляции сигналов с 1 и 2 датчика: ' ...
    num2str(max(R12_corr))]);
disp(['Максимум корреляции сигналов с 3 и 4 датчика: ' ...
    num2str(max(R34_corr))]);

clear f1 f2 R12 R34 f1_corr f2_corr R12_corr R34_corr shift