clear; close all; clc
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Номер дефекта

% Получаем путь к текущей папке
current_path = pwd;
% Разбиваем путь на отдельные части
[parent_path, ~] = fileparts(current_path);
% Получаем название предыдущей папки
[parent_path, lot_number] = fileparts(parent_path);
[~, road_name] = fileparts(parent_path);
% Путь к папке main
main_path = fileparts(mfilename('fullpath'));

clear parent_path current_path
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Загрузка данных

X1 = readtable('1.txt', 'HeaderLines', 1);
X2 = readtable('2.txt', 'HeaderLines', 1);
X3 = readtable('3.txt', 'HeaderLines', 1);
X4 = readtable('4.txt', 'HeaderLines', 1);

% Ускорения акселерометров на разных буксах в условных единицах
accel_1 = table2array(X1(:,4));
accel_2 = table2array(X2(:,4));
accel_3 = table2array(X3(:,4));
accel_4 = table2array(X4(:,4));

% Перевод в ускорения с помощью масштабного коэффицинта и сдвига нуля
koef = load(fullfile(main_path, 'p1_K_ADXL1001_Z.txt')); %%#ok<LOAD>
accel_1 = (accel_1 - koef(1, 1)) / koef(1, 2);
accel_2 = (accel_2 - koef(2, 1)) / koef(2, 2);
accel_3 = (accel_3 - koef(3, 1)) / koef(3, 2);
accel_4 = (accel_4 - koef(4, 1)) / koef(4, 2);

load('..\acc_1'); load('..\acc_2')
timeStamp_1 = acc_1(:,1); timeStamp_2 = acc_2(:,1);
sysCoord_1 =  acc_1(:,3); sysCoord_2 =  acc_2(:,3);

% Скорость на участке
V1_mean = mean((sysCoord_1(2:end) - sysCoord_1(1:end - 1)) .\ ...
    (timeStamp_1(2:end) - timeStamp_1(1:end - 1)));
V2_mean = mean((sysCoord_2(2:end) - sysCoord_2(1:end - 1)) .\ ...
    (timeStamp_2(2:end) - timeStamp_2(1:end - 1)));

clear ADXL1002_Z acc_1 acc_2 timeStamp_1 timeStamp_2 sysCoord_1 sysCoord_2
clear X1 X2 X3 X4 Fs
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Обработка сырых данных

% Приведение показаний акселерометрой в нулевой линии
accel_1 = accel_1 - mean(accel_1);
accel_2 = accel_2 - mean(accel_2);
accel_3 = accel_3 - mean(accel_3);
accel_4 = accel_4 - mean(accel_4);

% Нормирование ускорений на квадрат скорости
accel_1 = accel_1 / (V1_mean)^2;
accel_2 = accel_2 / (V2_mean)^2;
accel_3 = accel_3 / (V1_mean)^2;
accel_4 = accel_4 / (V2_mean)^2;

% accel_1 = accel_1 + 4 * randn(size(accel_1));
% accel_2 = accel_2 + 4 * randn(size(accel_2));
% accel_3 = accel_3 + 4 * randn(size(accel_3));
% accel_4 = accel_4 + 4 * randn(size(accel_4));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Общий графики участка ускорений акселерометров

figure(1)
plot(1:length(accel_1), accel_1,'k', 1:length(accel_2), accel_2, 'r', ...
    1:length(accel_3), accel_3, 'b', 1:length(accel_4), accel_4,'g')
    xlabel('Отсчеты')
    ylabel('Нормированные ускорения, м^{-1}')
%     ylabel('Ускорения, м/с^2')
    legend('1 Следом идущее','2 Впереди идущее', ...
        '3 Следом идущее','4 Впереди идущее')
%     legend('1 Впереди идущее L','2 Следом идущее L',
% '3 Впереди идущее R','4 Следом идущее R')
    grid on
    ylim([-20 28])
    xlim([0 length(accel_1)])
%     ylim([min(min(accel_1) - 6, min(accel_3) - 6) ...
%         max(max(accel_1) + 6,max(accel_3) + 6)])
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Корреляция сигналов до учета перераспределения энергии

% Огибающие сигнала для удобного вычисления максимума
[accel_1_up,~] = envelope(abs(accel_1), 100, 'peak');
[accel_2_up,~] = envelope(abs(accel_2), 100, 'peak');
[accel_3_up,~] = envelope(abs(accel_3), 100, 'peak');
[accel_4_up,~] = envelope(abs(accel_4), 100, 'peak');

% Нахождение индексов максимумов участков с повышенной амплитудой
if max(accel_1_up) >= max(accel_3_up)
    [~,index_1] = max(accel_1_up);
    [~,index_2] = max(accel_2_up);
    index_3 = index_1;
    index_4 = index_2;

elseif max(accel_1_up) <= max(accel_3_up)
    [~,index_3] = max(accel_3_up);
    [~,index_4] = max(accel_4_up);
    index_1 = index_3;
    index_2 = index_4;

end

clear accel_1_up accel_2_up accel_3_up accel_4_up

% Выбор длины участка
n = 1000; k = n + 499;

% Cовмещенные участки по максимумам
accel_1_comb = accel_1(index_1 - n:index_1 + k);
accel_2_comb = accel_2(index_2 - n:index_2 + k);
accel_3_comb = accel_3(index_1 - n:index_1 + k);
accel_4_comb = accel_4(index_2 - n:index_2 + k);

% Внесение поправок в индексы участков
[c_ind2, c_ind4] = index_correction(accel_1_comb, accel_2_comb, ...
    accel_3_comb, accel_4_comb);

% Индексы с учетом поправки
index_2 = index_2 + c_ind2;
index_4 = index_4 + c_ind4;

% Сигналы с учетом поправок
accel_1_comb = accel_1(index_1 - n:index_1 + k);
accel_2_comb = accel_2(index_2 - n:index_2 + k);
accel_3_comb = accel_3(index_3 - n:index_3 + k);
accel_4_comb = accel_4(index_4 - n:index_4 + k);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Графики ускорений и корреляций

% Графики совмещенных ускорений
figure(2)
plot(1:length(accel_1_comb), accel_1_comb,'k', ...
    1:length(accel_2_comb), accel_2_comb, 'r')
    xlabel('Отсчеты')
    ylabel('Нормированные ускорения, м^{-1}')
%     ylabel('Ускорения, м/с^2')
    legend('1 Следом идущее','2 Впереди идущее')
    grid on
    set(gca, 'YTick', -40:5:40)
    set(gca, 'XTick', 0:250:length(accel_3_comb))
    ylim([min(min(accel_1_comb) - 2, min(accel_3_comb) - 2) ...
        max(max(accel_1_comb) + 2, max(accel_3_comb) + 2)])

figure(3)
plot(1:length(accel_3_comb), accel_3_comb, 'b', ...
    1:length(accel_4_comb), accel_4_comb,'g')
    xlabel('Отсчеты')
    ylabel('Нормированные ускорения, м^{-1}')
%     ylabel('Ускорения, м/с^2')
    legend('3 Следом идущее','4 Впереди идущее')
    grid on  
    set(gca, 'YTick', -40:5:40)
    set(gca, 'XTick', 0:250:length(accel_3_comb))
    ylim([min(min(accel_1_comb) - 2, min(accel_3_comb) - 2) ...
        max(max(accel_1_comb) + 2, max(accel_3_comb) + 2)])

% Графики корреляции
figure(4)
s = 300;
[R_12,f1]=xcorr(accel_1_comb, accel_2_comb, 1*s, 'normalized');
% subplot(2,1,1)
plot(f1, R_12, 'b')
%     title("До учета")
    xlabel('Отсчеты')
    ylabel('Коэффициент корреляции')
    legend('1 и 2 колесо')
    ylim([-0.4, 1])
    xlim([-150 150])
    grid on
    set(gca, 'XTick', -300:25:300)
    set(gca, 'YTick',-0.6:0.2:1)
    
% s = 300;
[R_34,f2]=xcorr(accel_3_comb, accel_4_comb, 1*s, 'normalized');
% subplot(2,1,2)
% plot(f2, R_34)
%     xlabel('Отсчеты')
%     ylabel('Коэффициент корреляции')
%     legend('3 и 4 колесо')
%     ylim([-0.6, 1])
%     xlim([-150 150])
%     grid on
%     set(gca, 'XTick', -300:25:300)
%     set(gca, 'YTick',-0.5:0.2:1)

disp(['Максимум корреляции сигналов с 1 и 2 ' ...
    'датчика до учета: ' num2str(max(R_12))]);
disp(['Максимум корреляции сигналов с 3 и 4 ' ...
    'датчика до учета: ' num2str(max(R_34))]);

clear c_ind2 c_ind4 f1 f2 s R_12 R_34
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Учет коэффициентов перераспределения сигналов
K_axis_m = 0.27; K_rail_m = 0.21; K_diag_m = 0.20;
% Hd = filter_cheb;
Hd = filter_lowpass;
order = 334;
z = zeros(order / 2, 1);
% 
% accel_1_comb_filt = circshift(filter(Hd, [z; accel_1_comb; z]), ...
%     [-order/2, 0]);
% accel_2_comb_filt = circshift(filter(Hd, [z; accel_2_comb; z]), ...
%     [-order/2, 0]);
% accel_3_comb_filt = circshift(filter(Hd, [z; accel_3_comb; z]), ...
%     [-order/2, 0]);
% accel_4_comb_filt = circshift(filter(Hd, [z; accel_4_comb; z]), ...
%     [-order/2, 0]);
% 
% % Убираем отступы из сигналов
% accel_1_comb_filt = accel_1_comb_filt(order/2+1:end-order/2);
% accel_2_comb_filt = accel_2_comb_filt(order/2+1:end-order/2);
% accel_3_comb_filt = accel_3_comb_filt(order/2+1:end-order/2);
% accel_4_comb_filt = accel_4_comb_filt(order/2+1:end-order/2);

% a1 = accel_1_comb - accel_3_comb_filt * K_axis_m - ...
%     accel_2_comb_filt * K_rail_m - accel_4_comb_filt * K_diag_m;
% a2 = accel_2_comb - accel_4_comb_filt * K_axis_m - ... 
%     accel_1_comb_filt * K_rail_m - accel_3_comb_filt * K_diag_m;
% a3 = accel_3_comb - accel_1_comb_filt * K_axis_m - ...
%     accel_4_comb_filt * K_rail_m - accel_2_comb_filt * K_diag_m;
% a4 = accel_4_comb - accel_2_comb_filt * K_axis_m - ...
%     accel_3_comb_filt * K_rail_m - accel_1_comb_filt * K_diag_m;

% a1 = accel_1_comb - accel_3_comb_filt * K_axis_m;
% a2 = accel_2_comb - accel_4_comb_filt * K_axis_m;
% a3 = accel_3_comb - accel_1_comb_filt * K_axis_m;
% a4 = accel_4_comb - accel_2_comb_filt * K_axis_m;

a1 = accel_1_comb - accel_3_comb * K_axis_m;
a2 = accel_2_comb - accel_4_comb * K_axis_m;
a3 = accel_3_comb - accel_1_comb * K_axis_m;
a4 = accel_4_comb - accel_2_comb * K_axis_m;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Графики

figure(5)
plot(1:length(a1), a1,'k', 1:length(a2), a2, 'r')
    xlabel('Отсчеты')
    ylabel('Нормированные ускорения, м^{-1}')
%     ylabel('Ускорения, м/с^2')
    legend('1 Следом идущее','2 Впереди идущее')
    grid on  
    set(gca, 'YTick', -40:5:40)
    set(gca, 'XTick', 0:250:length(accel_3_comb))
    ylim([min(min(accel_1_comb) - 2, min(accel_3_comb) - 2) ...
        max(max(accel_1_comb) + 2, max(accel_3_comb) + 2)])

figure(6)
plot(1:length(a3), a3, 'b', 1:length(a4), a4,'g')
    xlabel('Отсчеты')
    ylabel('Нормированные ускорения, м^{-1}')
%     ylabel('Ускорения, м/с^2')
    legend('3 Следом идущее','4 Впереди идущее')
    grid on
    set(gca, 'YTick', -40:5:40)
    set(gca, 'XTick', 0:250:length(accel_3_comb))
    ylim([min(min(accel_1_comb) - 2, min(accel_3_comb) - 2) ...
        max(max(accel_1_comb) + 2, max(accel_3_comb) + 2)])

figure(7)
s = 300;
[R_12,f] = xcorr(a1, a2, 1*s, 'normalized');
% subplot(2,1,1)
plot(f, R_12, 'b')
%     title("После учета")
    xlabel('Отсчеты')
    ylabel('Коэффициент корреляции')
    legend('1 и 2 колесо')
    ylim([-0.4, 1])
    xlim([-150 150])
    grid on
    set(gca, 'XTick', -300:25:300)
    set(gca, 'YTick',-0.6:0.2:1)
    
[R_34,f] = xcorr(a3, a4, 1*s, 'normalized');
% subplot(2,1,2)
% plot(f, R_34)
%     xlabel('Отсчеты')
%     ylabel('Коэффициент корреляции')
%     legend('3 и 4 колесо')
%     ylim([-0.6, 1])
%     xlim([-150 150])
%     grid on 
%     set(gca, 'XTick', -300:25:300)
%     set(gca, 'YTick',-0.5:0.2:1)

disp(['Максимум корреляции сигналов с 1 и 2 ' ...
    'датчика после учета: ' num2str(max(R_12))]);
disp(['Максимум корреляции сигналов с 3 и 4 ' ...
    'датчика после учета: ' num2str(max(R_34))]);

clear f s R_12 R_34

%%
folder_path = pwd();
main_path = fileparts(mfilename('fullpath'));

% Загурзка данных ручных измерений
Manual_measure = load(fullfile(main_path, "manual_measurement.txt"));
[name_area, ~, ~] = fileparts(folder_path);
[~, idx_area, ~] = fileparts(name_area);
idx_area = str2double(idx_area) + 1;

accel_1_comb = accel_1_comb * (V1_mean)^2;
accel_2_comb = accel_2_comb * (V2_mean)^2;

a1 = a1 * (V1_mean)^2;
a2 = a2 * (V2_mean)^2;

Ts = 3.200005908651970e-05;
m = 650;
n = 1200;

% Глубина дефекта исходная в мм
h1 = (10^3) * cumsum(cumsum(accel_1_comb(m:n) * Ts) * Ts);
h2 = (10^3) * cumsum(cumsum(accel_1_comb(m:n) * Ts) * Ts);
hcorr1 = (10^3) * cumsum(cumsum(a1(m:n) * Ts) * Ts);
hcorr2 = (10^3) * cumsum(cumsum(a2(m:n) * Ts) * Ts);

figure
plot(h1, 'k')
hold on 
plot(hcorr1, 'b')
xlims = xlim;
line([xlims(1) xlims(2)], [-mean(Manual_measure(idx_area, :)) -mean(Manual_measure(idx_area, :))], 'Color', 'r'); 
hold off 
grid on
set(gca, 'XTick', -300:50:1000)
set(gca, 'YTick',-5:0.5:5)
ylabel('Глубина дефекта, мм')
xlabel('Отчеты')
legend('Исходные данные', 'Скорректированные данные', 'Ручное измерение')
% title('Колесо 1')
ylim([-2 1])

figure
plot(h2, 'k')
hold on 
plot(hcorr2, 'b')
xlims = xlim;
line([xlims(1) xlims(2)], [-mean(Manual_measure(idx_area, :)) -mean(Manual_measure(idx_area, :))], 'Color', 'r'); 
hold off
grid on
set(gca, 'XTick', -300:50:1000)
set(gca, 'YTick',-5:0.5:5)
ylabel('Глубина дефекта, мм')
% xlabel('Длинна дефекта, отчеты')
xlabel('Отчеты')
legend('Исходные данные', 'Скорректированные данные', 'Ручное измерение')
% title('Колесо 2')
ylim([-2 1])


%%
accel_1_comb_filt = accel_1_comb / (V1_mean)^2;
accel_2_comb_filt = accel_2_comb / (V2_mean)^2;
accel_3_comb_filt = accel_3_comb / (V1_mean)^2;
accel_4_comb_filt = accel_4_comb / (V2_mean)^2;

accel_1_comb = accel_1_comb / (V1_mean)^2;
accel_2_comb = accel_2_comb / (V2_mean)^2;
accel_3_comb = accel_3_comb / (V1_mean)^2;
accel_4_comb = accel_4_comb / (V2_mean)^2;


accel_1_comb_filt = circshift(filter(Hd, [z; accel_1_comb_filt; z]), ...
    [-order/2, 0]);
accel_2_comb_filt = circshift(filter(Hd, [z; accel_2_comb_filt; z]), ...
    [-order/2, 0]);
accel_3_comb_filt = circshift(filter(Hd, [z; accel_3_comb_filt; z]), ...
    [-order/2, 0]);
accel_4_comb_filt = circshift(filter(Hd, [z; accel_4_comb_filt; z]), ...
    [-order/2, 0]);

% Убираем отступы из сигналов
accel_1_comb_filt = accel_1_comb_filt(order/2+1:end-order/2);
accel_2_comb_filt = accel_2_comb_filt(order/2+1:end-order/2);
accel_3_comb_filt = accel_3_comb_filt(order/2+1:end-order/2);
accel_4_comb_filt = accel_4_comb_filt(order/2+1:end-order/2);


% accel_1_comb_filt = filter(Hd, accel_1_comb_filt);


% Графики совмещенных ускорений
figure
plot(1:length(accel_1_comb), accel_1_comb,'k', 1:length(accel_1_comb_filt), accel_1_comb_filt,'c')
    xlabel('Отсчеты')
    ylabel('Нормированные ускорения, м^{-1}')
%     ylabel('Ускорения, м/с^2')
    legend('Исходный','Фильтрованный')
    grid on
%     set(gca, 'YTick', -40:5:40)
    set(gca, 'XTick', 0:250:length(accel_3_comb))
%     ylim([min(min(accel_1_comb) - 2, min(accel_3_comb) - 2) ...
%         max(max(accel_1_comb) + 2, max(accel_3_comb) + 2)])

figure
s = 300;
[R_12,f] = xcorr(accel_1_comb, accel_1_comb_filt, 1*s, 'normalized');
% subplot(2,1,1)
plot(f, R_12, 'b')
%     title("После учета")
    xlabel('Отсчеты')
    ylabel('Коэффициент корреляции')
    legend('1 и 2 колесо')
%     ylim([-0.4, 1])
%     xlim([-150 150])
    grid on
    set(gca, 'XTick', -300:25:300)
    set(gca, 'YTick',-0.6:0.2:1)