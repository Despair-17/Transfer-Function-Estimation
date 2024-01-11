function [signal_input, signal_output_axis, ...
                signal_output_rail, signal_output_diag, Ts] = ...
                processAndCombine(folder_path, main_path, len_before, ...
                len_after, wheel_number)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Детальное описание функции:
%   processAndCombine - принимает на вход конкретную директорию из которой 
%   необходимо загрузить данные и wheel_number, после их загрузки с помощью
%   функции loadAndPreprocess объединяет их в массивы всех возможных вари-
%   антов входных и выходных сигналов системы, а после возвращает их, длина
%   возвращаемых массивов равна длине одного исследуемого участка умножен-
%   ное на 4 
%
% Входные аргументы:
% - folder_path - путь к папке c сырыми файлами показаний акселерометров, 
%   которые установлены на буксах тележки вагона.
% - main_path - путь к папке main, где находятся остальные объекты (другие
%   используемые скрипты, функции и фильтр).
% - len_before, len_after - эти аргументы используются для определения
%   длинн исследуемых участков.
%
% Выходные аргументы:
% - signal_input - последовательно объединенные все входные сигналы
%   системы.
% - signal_output_axis - последовательно объединенные все выходные сигналы
%   системы (колесо на которое воздействовали -> колесо той же колесной 
%   пары).
% - signal_output_rail - последовательно объединенные все выходные сигналы
%   системы (колесо на которое воздействовали -> колесу находящемся на том 
%   же рельсе).
% - signal_output_diag - последовательно объединенные все выходные сигналы
%   системы (колесо на которое воздействовали -> колесо находящиеся по 
%   диагонали).
% - Ts - средний период дискретизации для четырех файлов с данными.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Код функции начинается здесь

    % Вызов функции которая загружает данные и производит обработку
    [accel, ~, ~, Ts] = ...
        loadAndPreprocess(folder_path, main_path, wheel_number);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Выделение участков с полезными сигналами
    
    % Огибающие сигнала для удобного вычисления максимума
    [~, idx] = max(envelope(abs(accel(:, wheel_number)), 100, 'peak'));

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Фильтрация сигналов

    [Hd, order] = filter_lowpass;
    z = zeros(order / 2, 1);

    % Наполнение нулями, чтобы сохранить положение сигнала во времени, а 
    % потом его совмещение с изначальным положением

    % Сдвигаем выходной сигнал на половину порядка фильтра
    accel_filt = [...
        circshift(filter(Hd, [z; accel(:, 1); z]), [-order/2, 0]), ...
        circshift(filter(Hd, [z; accel(:, 2); z]), [-order/2, 0]), ...
        circshift(filter(Hd, [z; accel(:, 3); z]), [-order/2, 0]), ...
        circshift(filter(Hd, [z; accel(:, 4); z]), [-order/2, 0])  ...
                 ];

    % Убираем отступы из сигналов
    accel_filt = accel_filt(order/2+1:end-order/2, :);

    clear z Hd order 
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Организация данных

    % Выбор длин исследуемых участков
    z = zeros(len_before + len_after + 1, 1);
    len_signsl = idx - len_before:idx + len_after;

    if wheel_number == 1
        % Когда воздействие происходит на колеса 1

        % Входные и выходнные сигналы дополненные нулями
        signal_input = [z; accel_filt(len_signsl, 1)];
        signal_output_axis = [z; accel_filt(len_signsl, 3)];
        signal_output_rail = [z; accel_filt(len_signsl, 2)];
        signal_output_diag = [z; accel_filt(len_signsl, 4)];

    elseif wheel_number == 2

        % Когда воздействие происходит на колеса 2
        signal_input = [z; accel_filt(len_signsl, 2)];
        signal_output_axis = [z; accel_filt(len_signsl, 4)];
        signal_output_rail = [z; accel_filt(len_signsl, 1)];
        signal_output_diag = [z; accel_filt(len_signsl, 3)];
   
    elseif wheel_number == 3

        % Когда воздействие происходит на колеса 3
        signal_input = [z; accel_filt(len_signsl, 3)];
        signal_output_axis = [z; accel_filt(len_signsl, 1)];
        signal_output_rail = [z; accel_filt(len_signsl, 4)];
        signal_output_diag = [z; accel_filt(len_signsl, 2)];

    elseif wheel_number == 4

        % Когда воздействие происходит на колеса 4
        signal_input = [z; accel_filt(len_signsl, 4)];
        signal_output_axis = [z; accel_filt(len_signsl, 2)];
        signal_output_rail = [z; accel_filt(len_signsl, 3)];
        signal_output_diag = [z; accel_filt(len_signsl, 1)];

    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Код функции заканчивается здесь
end


