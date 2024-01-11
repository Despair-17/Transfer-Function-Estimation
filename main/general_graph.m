clear all
close all
clc
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Номер дефекта

% Папка куда сохранять гарфики
folderName = 'D:\programs\common\plots\1002';

% Получаем путь к текущей папке
current_path = pwd;
% Разбиваем путь на отдельные части
[parent_path, ~] = fileparts(current_path);
% Получаем название предыдущей папки
[parent_path, lot_number] = fileparts(parent_path);
[~, road_name] = fileparts(parent_path);
% Получаем путь к основному скрипту
main_path = fileparts(mfilename('fullpath'));

disp("Проезд " + road_name + " участок " + lot_number)

% clear parent_path current_path
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Путь
for i = 0:33
path1 = string(current_path) + "\" + num2str(i) + "\raw\";
path2 = string(current_path) + "\" + num2str(i) + "\";

% Загрузка данных

% load('ADXL1002_Z.mat')
X1 = readtable(path1 + '1.txt', 'HeaderLines', 1);
X2 = readtable(path1 + '2.txt', 'HeaderLines', 1);
X3 = readtable(path1 + '3.txt', 'HeaderLines', 1);
X4 = readtable(path1 + '4.txt', 'HeaderLines', 1);

% Ускорения акселерометров на разных буксах в условных единицах
accel_1 = table2array(X1(:, 7));
accel_2 = table2array(X2(:, 7));
accel_3 = table2array(X3(:, 7));
accel_4 = table2array(X4(:, 7));

% Перевод в ускорения с помощью масштабного коэффицинта и сдвига нуля
koef = load(fullfile(main_path, 'p1_K_ADXL1002_Z.txt')); %%#ok<LOAD>

accel_1 = (accel_1 - koef(1, 1)) / koef(1, 2); 
accel_2 = (accel_2 - koef(2, 1)) / koef(2, 2);
accel_3 = (accel_3 - koef(3, 1)) / koef(3, 2);
accel_4 = (accel_4 - koef(4, 1)) / koef(4, 2);

n = 8;
accel_1 = smoothdata(accel_1, 'movmean', n);
accel_2 = smoothdata(accel_2, 'movmean', n);
accel_3 = smoothdata(accel_3, 'movmean', n);
accel_4 = smoothdata(accel_4, 'movmean', n);

load(path2 + '\acc_1')
load(path2 + '\acc_2')

timeStamp_1 = acc_1(:,1);
sysCoord_1 =  acc_1(:,3);
timeStamp_2 = acc_2(:,1);
sysCoord_2 =  acc_2(:,3);

% Скорость на участке
V1_mean = mean((sysCoord_1(2:end) - sysCoord_1(1:end - 1)) .\ ...
    (timeStamp_1(2:end) - timeStamp_1(1:end - 1)));
V2_mean = mean((sysCoord_2(2:end) - sysCoord_2(1:end - 1)) .\ ...
    (timeStamp_2(2:end) - timeStamp_2(1:end - 1)));

% Частота дискретизации
% Fs = length(acc_1(:,5)) / abs(((timeStamp_2(end) /10^8) - ...
%     (timeStamp_1(1) /10^8)));

clear ADXL1002_Z acc_1 acc_2 timeStamp_1 timeStamp_2 sysCoord_1 sysCoord_2
clear X1 X2 X3 X4 Fs
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Обработка сырых данных

% Приведение показаний акселерометрой в нулевой линии
accel_1 = accel_1 - mean(accel_1);
accel_2 = accel_2 - mean(accel_2);
accel_3 = accel_3 - mean(accel_3);
accel_4 = accel_4 - mean(accel_4);

% Нормирование ускорений на квадрат скорости
% accel_1 = accel_1 / (V1_mean)^2;
% accel_2 = accel_2 / (V2_mean)^2;
% accel_3 = accel_3 / (V1_mean)^2;
% accel_4 = accel_4 / (V2_mean)^2;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Общий графики участка ускорений акселерометров

figure(1)
plot(1:length(accel_1), accel_1,'k', 1:length(accel_2), accel_2, 'r', ...
    1:length(accel_3), accel_3, 'b', 1:length(accel_4), accel_4,'g')
    xlabel('Отсчеты')
%     ylabel('Нормированные ускорения, м^-1')
    ylabel('Ускорения, м/с^2')
    legend('1 Следом идущее','2 Впереди идущее', ...
        '3 Следом идущее','4 Впереди идущее')
%     legend('1 Впереди идущее','2 Следом идущее', ...
%           '3 Впереди идущее','4 Следом идущее')
    grid on
    xlim([0 length(accel_1)])
    ylim([min(min(accel_1) - 2, min(accel_3) - 2) ...
        max(max(accel_1) + 2, max(accel_3) + 2)])

if ~isfolder(folderName)
    mkdir(folderName);
end

filename = string(num2str(i)) + '.png';
fullFileName = fullfile(folderName, filename);
saveas(gcf, fullFileName, 'png');

end