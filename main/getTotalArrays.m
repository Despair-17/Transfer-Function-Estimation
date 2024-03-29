function [data_input_outputs_signals, Ts] ...
          = getTotalArrays(folder_travel, folder_area, wheel_number)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Детальное описание функции:
%   getTotalArrays - принимает на вход название папок проездов и названия
%   определенных участков и возвращает объеденные массив данных различных
%   входных сигналов системы и всеми возможными выходными данными, размер
%   возаращаемого массива зависит от количества переданых названий папок, 
%   где содержатся сырые данные.
%
% Входные аргументы:
% - folder_travel - название папок проездов.
% - folder_area - название папок определенных участков.
%
% Выходные аргументы:
% - data_input_outputs_signals - общий массив всех исследуемых данных.
%   Первый столбец - объединенный массив данных входов системы. Второй
%   столбец - последовательно объединенные все выходные сигналы системы 
%   (колесо на которое воздействовали -> колесо той же колесной пары).
%   Третий столбец - последовательно объединенные все выходные сигналы 
%   системы (колесо на которое воздействовали -> колесу находящемся 
%   на том же рельсе). Четвертый столбец - последовательно объединенные 
%   все выходные сигналы системы (колесо на которое воздействовали -> 
%   колесо находящиеся по диагонали).
% - Ts - двухуровневое среднее значение периода дискретизации.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Код функции начинается здесь

    % Получаем путь к основному скрипту
    main_path = fileparts(mfilename('fullpath'));

    % Формируем путь к папке с данными
    data_path = fullfile(main_path, '..', 'data');

    % Задают длинную исследуемого участка
    len_before = 500; 
    len_after = 1499;
  
    % Количество обрабатываемых объектов
    num_sections = sum(cellfun(@numel, folder_area));

    % Создание пустых массивов нужной длинны для всех данных
    len_chunk = (len_before + len_after + 1) * 2;
    data_input_outputs_signals = zeros(len_chunk * num_sections, 4);
    Ts = zeros(num_sections, 1);
    
    % Индекс для записи данных в массив
    idx = 1;
    idx_Ts = 1;
    
    fprintf('Запуск обхода каталогов, в сумме их %d штук.\n', ...
        num_sections);

    % Внешний цикл проходится по индекса массива с папками folder_travel
    for i = 1:numel(folder_travel)

        % Путь к проездам (folder_travel)
        travel_path = fullfile(data_path, folder_travel{i});
        fprintf('Обход каталогов внутри: %s...\n', folder_travel{i});

        for j = 1:numel(folder_area{i})
            
            % Путь к сырым данным
            raw_path = fullfile(travel_path, folder_area{i}{j}, 'raw');

            % Запуск функции для сбора данных с определенного участка
            [
             data_input_outputs_signals(idx:idx + len_chunk - 1, 1), ...
             data_input_outputs_signals(idx:idx + len_chunk - 1, 2), ...
             data_input_outputs_signals(idx:idx + len_chunk - 1, 3), ...
             data_input_outputs_signals(idx:idx + len_chunk - 1, 4), ...
             Ts(idx_Ts, 1) ...
            ]...
             = processAndCombine(raw_path, main_path, ...
                                 len_before, len_after, wheel_number);

            idx = idx + len_chunk;
            idx_Ts = idx_Ts + 1;

        end

    end
    
    % Среднее значение периода дискритизации
    Ts = mean(Ts);

    fprintf('Данные системы успешно записанны в оперативную память.\n');
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Код функции заканчивается здесь
end












