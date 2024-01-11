function [accel, V1mean, V2mean, Ts] = ...
                    loadAndPreprocess(folder_path, main_path, wheel_number)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Детальное описание функции:
%   loadAndPreprocess - функция загружает сырые данные с указанной папки и
%   возвращет массив данных после обработки. Обработка включает перевод ус-
%   ловных единиц данных в ускорения, приведения показаний к нулевой линии,
%   а также нормирование ускорений по скорости движения.
%
% Входные аргументы:
%   folder_path - путь к сырым данным.
%   main_path - путь к папке main.
%
% Выходные аргументы:
%   accel - массив данных нормированных ускорений.
%   V1mean, V2mean - скорости движения в момент воздействия
%   Ts - период дискретизации.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Код функции начинается здесь

    % - Загрузка сырых данных с txt файлов.
    X1 = readtable(fullfile(folder_path, '1.txt'), 'HeaderLines', 1);
    X2 = readtable(fullfile(folder_path, '2.txt'), 'HeaderLines', 1);
    X3 = readtable(fullfile(folder_path, '3.txt'), 'HeaderLines', 1);
    X4 = readtable(fullfile(folder_path, '4.txt'), 'HeaderLines', 1);

    % Ускорения акселерометров в условных единицах
    min_len = min([height(X1), height(X2), height(X3), height(X4)]);
    accel = [table2array(X1(1:min_len, 4)), ...
             table2array(X2(1:min_len, 4)), ...
             table2array(X3(1:min_len, 4)), ...
             table2array(X4(1:min_len, 4))];

    accel2 = [table2array(X1(1:min_len, 9)), ...
             table2array(X2(1:min_len, 9)), ...
             table2array(X3(1:min_len, 9)), ...
             table2array(X4(1:min_len, 9))];

    % Средний период дискретизации
    Ts = mean([mean(diff(table2array(X1(1:min_len, 1)) / 10^8)), ...
                mean(diff(table2array(X2(1:min_len, 1)) / 10^8)), ...
                mean(diff(table2array(X3(1:min_len, 1)) / 10^8)), ...
                mean(diff(table2array(X4(1:min_len, 1)) / 10^8))]);

    clear X1 X2 X3 X4

    % Перевод в ускорения с помощью масштабного коэффицинта и сдвига нуля
    koef = load(fullfile(main_path, 'p1_K_ADXL1001_Z.txt')); %%#ok<LOAD>
    accel(:, 1) = (accel(:, 1) - koef(1, 1)) / koef(1, 2);
    accel(:, 2) = (accel(:, 2) - koef(2, 1)) / koef(2, 2);
    accel(:, 3) = (accel(:, 3) - koef(3, 1)) / koef(3, 2);
    accel(:, 4) = (accel(:, 4) - koef(4, 1)) / koef(4, 2);

    % Перевод в ускорения с помощью масштабного коэффицинта и сдвига нуля
    koef = load(fullfile(main_path, 'p1_K_ADXL354_X.txt')); %%#ok<LOAD>
    accel2(:, 1) = (accel2(:, 1) - koef(1, 1)) / koef(1, 2);
    accel2(:, 2) = (accel2(:, 2) - koef(2, 1)) / koef(2, 2);
    accel2(:, 3) = (accel2(:, 3) - koef(3, 1)) / koef(3, 2);
    accel2(:, 4) = (accel2(:, 4) - koef(4, 1)) / koef(4, 2);

    % Замена нужных ускорений на более чувствительный датчик
%     if wheel_number == 1
% 
%         accel = [accel(:, 1), accel2(:, 2:4)];
% 
%     elseif wheel_number == 2
% 
%         accel = [accel2(:, 1), accel(:, 2), accel2(:, 3:4)];
% 
%     elseif wheel_number == 3
% 
%         accel = [accel2(:, 1:2), accel(:, 3), accel2(:, 4)];
% 
%     elseif wheel_number == 4
% 
%         accel = [accel2(:, 1:3), accel(:, 4)];
% 
%     end

    % Загрузка временных меток и систмных координат
    [parent_dir, ~, ~] = fileparts(folder_path);
    acc_1 = load(fullfile(parent_dir, 'acc_1'));
    acc_2 = load(fullfile(parent_dir, 'acc_2'));

    timeStamp_1 = acc_1(:,1);
    timeStamp_2 = acc_2(:,1);
    sysCoord_1 =  acc_1(:,3);
    sysCoord_2 =  acc_2(:,3);

    % Расчет средней скорости на участке
    V1mean = mean((sysCoord_1(2:end) - sysCoord_1(1:end - 1)) .\ ...
        (timeStamp_1(2:end) - timeStamp_1(1:end - 1)));
    V2mean = mean((sysCoord_2(2:end) - sysCoord_2(1:end - 1)) .\ ...
        (timeStamp_2(2:end) - timeStamp_2(1:end - 1)));

    clear acc_1 acc_2 timeStamp_1 timeStamp_2 sysCoord_1 sysCoord_2 koef
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% Обработка сигналов

    % Приведение показаний акселерометрой к нулевой линии
    accel(:, 1) = (accel(:, 1) - mean(accel(:, 1)));
    accel(:, 2) = (accel(:, 2) - mean(accel(:, 2)));
    accel(:, 3) = (accel(:, 3) - mean(accel(:, 3)));
    accel(:, 4) = (accel(:, 4) - mean(accel(:, 4)));

    % Нормирование ускорений акселерометров на квадрат скорости
    accel(:, 1) = accel(:, 1) / (V1mean)^2;
    accel(:, 2) = accel(:, 2) / (V2mean)^2;
    accel(:, 3) = accel(:, 3) / (V1mean)^2;
    accel(:, 4) = accel(:, 4) / (V2mean)^2;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Код функции заканчивается здесь
end