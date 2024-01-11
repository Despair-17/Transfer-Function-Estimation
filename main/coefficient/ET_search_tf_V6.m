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
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Общий графики участка ускорений акселерометров

figure(1)
plot(1:length(accel_1), accel_1,'k', 1:length(accel_2), accel_2, 'r', ...
    1:length(accel_3), accel_3, 'b', 1:length(accel_4), accel_4,'g')
    xlabel('Отсчеты')
    ylabel('Нормированные ускорения, м^-1')
%     ylabel('Ускорения, м/с^2')
    legend('1 Следом идущее L','2 Впереди идущее L', ...
        '3 Следом идущее R','4 Впереди идущее R')
%     legend('1 Впереди идущее L','2 Следом идущее L',
% '3 Впереди идущее R','4 Следом идущее R')
    grid on
    xlim([0 length(accel_1)])
    ylim([min(min(accel_1) - 5, min(accel_3) - 5) ...
        max(max(accel_1) + 5,max(accel_3) + 5)])
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
    ylabel('Нормированные ускорения, м^-1')
%     ylabel('Ускорения, м/с^2')
    legend('1 Следом идущее L','2 Впереди идущее L')
    grid on  
    ylim([min(min(accel_1_comb) - 1, min(accel_3_comb) - 1) ...
        max(max(accel_1_comb) + 1, max(accel_3_comb) + 1)])

figure(3)
plot(1:length(accel_3_comb), accel_3_comb, 'b', ...
    1:length(accel_4_comb), accel_4_comb,'g')
    xlabel('Отсчеты')
    ylabel('Нормированные ускорения, м^-1')
%     ylabel('Ускорения, м/с^2')
    legend('3 Следом идущее R','4 Впереди идущее R')
    grid on  
    ylim([min(min(accel_1_comb) - 1, min(accel_3_comb) - 1) ...
        max(max(accel_1_comb) + 1, max(accel_3_comb) + 1)])

% Графики корреляции
figure(4)
s = 300;
[R_12,f1]=xcorr(accel_1_comb, accel_2_comb, 1*s, 'normalized');
subplot(2,1,1)
plot(f1, R_12)
    title("До учета")
    xlabel('Отсчеты')
    ylabel('Коэффициент корреляции')
    legend('1 и 2 колесо')
    ylim([-0.7, 1])
    grid on
    set(gca, 'YTick',-0.5:0.2:1)
    
s = 300;
[R_34,f2]=xcorr(accel_3_comb, accel_4_comb, 1*s, 'normalized');
subplot(2,1,2)
plot(f2, R_34)
    xlabel('Отсчеты')
    ylabel('Коэффициент корреляции')
    legend('3 и 4 колесо')
    ylim([-0.7, 1])
    grid on
    set(gca, 'YTick',-0.5:0.2:1)

disp(['Максимум корреляции сигналов с 1 и 2 ' ...
    'датчика до учета: ' num2str(max(R_12))]);
disp(['Максимум корреляции сигналов с 3 и 4 ' ...
    'датчика до учета: ' num2str(max(R_34))]);

clear c_ind2 c_ind4 f1 f2 s R_12 R_34
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Фильтрация сигналов
% Hd = filter_main; % - для сравнения
Hd = filter_lowpass;
order = 334;
z = zeros(order / 2, 1);

% Наполнение нулями, чтобы сохранить положение сигнала во времени, а потом
% его совмещение с изначальным положением

% Сдвигаем выходной сигнал на половину порядка фильтра
accel_1_comb_filt = circshift(filter(Hd, [z; accel_1_comb; z]), ...
    [-order/2, 0]);
accel_2_comb_filt = circshift(filter(Hd, [z; accel_2_comb; z]), ...
    [-order/2, 0]);
accel_3_comb_filt = circshift(filter(Hd, [z; accel_3_comb; z]), ...
    [-order/2, 0]);
accel_4_comb_filt = circshift(filter(Hd, [z; accel_4_comb; z]), ...
    [-order/2, 0]);

% Убираем отступы из сигналов
accel_1_comb_filt = accel_1_comb_filt(order/2+1:end-order/2);
accel_2_comb_filt = accel_2_comb_filt(order/2+1:end-order/2);
accel_3_comb_filt = accel_3_comb_filt(order/2+1:end-order/2);
accel_4_comb_filt = accel_4_comb_filt(order/2+1:end-order/2);

figure(5)
plot(1:length(accel_1_comb_filt), accel_1_comb_filt,'k', ...
    1:length(accel_2_comb_filt), accel_2_comb_filt, 'r')
    xlabel('Отсчеты')
    ylabel('Нормированные ускорения, м^-1')
%     ylabel('Ускорения, м/с^2')
    legend('1 Следом идущее L','2 Впереди идущее L')
    grid on  
    ylim([min(min(accel_1_comb) - 1, min(accel_3_comb) - 1) ...
        max(max(accel_1_comb) + 1, max(accel_3_comb) + 1)])

figure(6)
plot(1:length(accel_3_comb_filt), accel_3_comb_filt, 'b', ...
    1:length(accel_4_comb_filt), accel_4_comb_filt,'g')
    xlabel('Отсчеты')
    ylabel('Нормированные ускорения, м^-1')
%     ylabel('Ускорения, м/с^2')
    legend('3 Следом идущее R','4 Впереди идущее R')
    grid on 
    ylim([min(min(accel_1_comb) - 1, min(accel_3_comb) - 1) ...
        max(max(accel_1_comb) + 1, max(accel_3_comb) + 1)])

clear z order Hd 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
%% Оценка передаточных функций системы

Fs = 31250;
Ts = 1/Fs;
t = 0:1/Fs:(length(accel_1_comb_filt) - 1) / Fs;

% создаем объект данных
data = iddata(accel_3_comb_filt, accel_1_comb_filt, Ts);
% оцениваем передаточную функцию второго порядка с одним выходом

prompt = 'Выберите вариант:';
option1 = 'Рассчитать передаточную функцию и сохранить';
option2 = 'Загрузить передаточную функцию с папки main';

% choice = menu(prompt, option1, option2);
% switch choice
%     case 1
%         choice_user = "y";
%     case 2
%         choice_user = "n";
% end
% Ts = 3.200005908651970e-05;
% if "y" == choice_user
    tf_axis = tfest(data, 8, 'Ts', Ts);
%     save("D:\Photo and Video and Other\Важно (учеба)\Maga\" + ...
%         "Диплом Магистра\semestr 3\Нарезка-Удомли\main\TF_axis.mat", "tf_axis")
% elseif "n" == choice_user
%     load("tf_axis.mat")
% end
% 
clear option1 option2 order prompt choice_user choice_user
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Итоговые результаты и графики
accel_3_restored = lsim(tf_axis, accel_1_comb_filt, t);

figure
plot(1:length(accel_3_comb_filt), accel_3_comb_filt,'k', 1:length(accel_3_restored), accel_3_restored,'r')
    xlabel('Отсчеты')
    ylabel('Нормированные ускорения, м^{-1}')
%     ylabel('Ускорения, м/с^2')
    legend('Исходный','Модель')
    grid on
    xlim([500 1750])
    ylim([-3 3])
%     set(gca, 'YTick', -40:5:40)
    set(gca, 'XTick', 0:250:length(accel_3_comb))
%     ylim([min(min(accel_1_comb) - 2, min(accel_3_comb) - 2) ...
%         max(max(accel_1_comb) + 2, max(accel_3_comb) + 2)])
% % отображение АЧХ и ФЧХ в децибелах
% figure(7)
% bode(tf_axis); 
%%

max_window_size = 32; % максимальный размер окна
max_overlap = 16; % максимальное перекрытие

sys_estimate = tfestimate(accel_3_comb_filt,accel_1_comb_filt, max_window_size, max_overlap, 2^14, Fs);
sys = tf(sys_estimate);

y_out = lsim(sys, accel_1_comb_filt, t);
y_out_column = y_out(:, end);

figure(10)
plot(accel_3_comb_filt, 'k')
hold on 
plot(y_out_column, 'r')
hold off
    xlabel('Отсчеты')
    ylabel('Нормированные ускорения, м^{-1}')
%     ylabel('Ускорения, м/с^2')
    legend('Исходный','Модель')
    grid on
    xlim([500 1750])
    ylim([-3 3])
%     set(gca, 'YTick', -40:5:40)
    set(gca, 'XTick', 0:250:length(accel_3_comb))

% % получаем выходной сигнал, зная входной
% if max(accel_1_comb_filt) > max(accel_3_comb_filt)
%     accel_1_restored = lsim(tf_axis, accel_1_comb_filt, t);
%     accel_1_restored = [accel_1_restored(55:end); zeros(55, 1)];
%     accel_2_restored = lsim(tf_axis, accel_2_comb_filt, t);
%     accel_2_restored = [accel_2_restored(55:end); zeros(55, 1)];
% 
% elseif  max(accel_1_comb_filt) < max(accel_3_comb_filt)
%     accel_1_restored = lsim(tf_axis, accel_3_comb_filt, t);
%     accel_1_restored = [accel_1_restored(55:end); zeros(55, 1)];
%     accel_2_restored = lsim(tf_axis, accel_4_comb_filt, t);
%     accel_2_restored = [accel_2_restored(55:end); zeros(55, 1)];
% 
% end
% %
% figure
% compare(data, tf_axis);
% title("axis")

% figure(8)
% % subplot(1,2,1)
% plot(1:length(accel_1_restored), accel_1_restored,'k', ...
%     1:length(accel_3_comb_filt), accel_3_comb_filt, 'r')
%     xlabel('Отсчеты')
%     ylabel('Нормированные ускорения, м^-1')
% %     ylabel('Ускорения, м/с^2')
% %     legend('1 Следом идущее L','2 Впереди идущее L')
%     legend('Восстановленный сигнал','Сигнал с датчика L')
%     grid on  
% 
% figure(9)
% % subplot(1,2,2)
% plot(1:length(accel_2_restored), accel_2_restored,'k', ...
%     1:length(accel_4_comb_filt), accel_4_comb_filt, 'r')
%     xlabel('Отсчеты')
%     ylabel('Нормированные ускорения, м^-1')
% %     ylabel('Ускорения, м/с^2')
% %     legend('1 Следом идущее L','2 Впереди идущее L')
%     legend('Восстановленный сигнал','Сигнал с датчика L')
%     grid on  

%     %% метод 2
% window_size = 1200; % размер окна
% overlap = 400; % перекрытие
% u1 = [accel_1_comb_filt; accel_2_comb_filt];
% y1 = [accel_3_comb_filt; accel_4_comb_filt];
% sys_estimate = tfestimate(accel_1_comb_filt, accel_3_comb_filt, window_size, overlap, 2500, Fs);
% sys = tf(sys_estimate);
% 
% t = linspace(0, 10, length(accel_1_comb_filt));
% y_out = lsim(sys, accel_1_comb_filt, t);
% y_out_column = y_out(:, end);
% 
% figure(100)
% plot(y_out_column)



% sys_estimate = tfestimate(accel_1_comb_filt, accel_3_comb_filt, 64, 32, 2500);
% [num, den] = tfdata(tf(sys_estimate));
% 
% sys = tf(num, den);
% 
% y = lsim(sys, accel_1_comb_filt, t);

% figure(100)
% plot(y)

% sys_estimate = tfestimate(accel_1_comb_filt, accel_3_comb_filt, window_size, overlap, window_size); % оценка передаточной функции
% num = real(sys_estimate)';
% den = real([1; -sys_estimate])';
% sys = tf(num, den);
% y = lsim(sys, accel_1_comb_filt, t); % выходной сигнал
