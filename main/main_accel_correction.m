clear; close all; clc; tic;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Детальное описание основного скрипта программы:
%   main_accel_correction - скрипт реализует основную логику учета перера-
%   спределения энерги, запаская функцию processFolders которая обходит
%   указанные проезды folder_travel и конкретные участки folder_area, в
%   свою очередь функция processFolders запускает функцию energyRedistri-
%   bution которая обрабатывает данные после чего сохраняет скорректирован-
%   ные ускорения в каталоги с исходными. В переменной wheel_choices нужно
%   указать перераспределения от каких колес нужно учесть (1 - Да, 0 - Нет)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Ввод данных для коррекции

folder_travel = {'нов-обратно', 'нов-туда', 'стар-обратно', 'стар-туда'};
folder_area = {{'0', '25', '26', '20', '21', '24', '27'}, ...
               {'0', '25', '26', '20', '21', '24', '27'}, ...
               {'0', '25', '26', '20', '21', '24', '27'}, ...
               {'0', '25', '26', '20', '21', '24', '27'}};

wheel_choices = {{[1 1 1 1], [1 1 1 1], [1 1 1 1], [1 1 1 0], ...
                  [1 1 1 0], [1 1 1 1], [1 1 1 1]}, ... 
                 {[1 1 1 1], [1 1 1 1], [1 1 1 1], [1 1 1 0], ...
                  [1 1 1 0], [1 1 1 1], [1 1 1 1]}, ...
                 {[1 1 1 1], [1 1 1 1], [1 1 1 1], [0 1 1 1], ...
                  [0 1 1 1], [0 1 1 1], [1 1 1 1]}, ...
                 {[1 1 1 1], [1 1 1 1], [1 1 1 1], [0 1 1 1], ...
                  [0 1 1 1], [1 1 1 1], [1 1 1 1]}};

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
flag = exist('folder_travel', 'var') & ...
       exist('folder_area', 'var') & ...
       exist('wheel_choices', 'var');

if ~flag
    disp(['Ошибка! Проверьте наличие в рабочем пространстве ' ...
        'переменных: folder_travel, folder_area, wheel_number.']);
    return
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fprintf('Запуск программы...\n')
fprintf(['Если скорректированные данные уже содержатся в указанных ' ...
    'папка они будет перезаписаны, либо можете прервать запуск.\n'])

prompt = 'Выбор варианта действий.';
option1 = 'Перезаписать скорректированные данные.';
option2 = 'Прервать выполнения программы.';

choice = menu(prompt, option1, option2);

switch choice
    case 1

    case 2
        fprintf('Запуск прерван.\n')
        return
end

clear option1 option2 prompt flag choice
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Путь к основному скрипту и данным для оценки передаточных функций
main_path = fileparts(mfilename('fullpath'));

% Создание пустой структуры с передаточными функциями
sys_struct = struct();

% Загрузка данных из файлов и сохранение в структуру
for i = 1:4
    filename = fullfile(main_path, ['../data/tf system/sysd_axis', ...
        num2str(i) ,'.mat']);
    sysd_axis = load(filename).sysd_axis;
    sys_struct.(['sysd_axis', num2str(i)]) = sysd_axis;

    filename = fullfile(main_path, ['../data/tf system/sysd_rail', ...
        num2str(i) ,'.mat']);
    sysd_rail = load(filename).sysd_rail;
    sys_struct.(['sysd_rail', num2str(i)]) = sysd_rail;

    filename = fullfile(main_path, ['../data/tf system/sysd_diag', ...
        num2str(i) ,'.mat']);
    sysd_diag = load(filename).sysd_diag;
    sys_struct.(['sysd_diag', num2str(i)]) = sysd_diag;
end

clear sysd_axis sysd_rail sysd_diag i
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Запуск обхода директорий

processFolders(folder_travel, folder_area, wheel_choices, sys_struct)