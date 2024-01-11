clear; close all; clc; tic;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Детальное описание основного скрипта программы:
%   main_tf_evaluation - скрипт организует основную логику оценивания пере-
%   даточных передаточных функций системы вагона тележки. Передаточная фун-
%   кция sysd_axis - это дискретная модель системы колеса на которое возде-
%   йствовали, колесо тойже колесной пары. sysd_rail - это дискретная моде-
%   ль системы колеса на которое воздействовали, колесо тойже колесной па-
%   ры. sysd_diag - это дискретная модель системы колеса на которое воздей-
%   ствовали, диагональное колесо. Скрипт вызывает функции getTotalArrays,
%   processAndCombine, evalSysTF. Первая собирает данные в один общий мас-
%   сив данных обходя по всем указанным папкам folder_travel, folder_area,
%   вторая собирает данные с одного участка, третья функция оценивает диск-
%   ретную модель системы для более подробного описания см. их описание.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Ввод данных для исследования

% Название папок проездов и определенных участков
% Директории для 1 и 2 колеса
folder_travel = {'нов-обратно', 'нов-туда'};
folder_area = {{'4', '6', '7', '8', '11', '17', '18', '19'}, ...
               {'4', '8', '9', '10', '11', '14', '15', '17', '19'}};
% wheel_number = 1;
wheel_number = 2;

% Директории для 3 и 4 колеса
% folder_travel = {'стар-обратно', 'стар-туда'};
% folder_area = {{'2', '3', '4', '5', '6', '7', '8', '9', '10', '11', ...
%                   '16', '17', '18', '19'}, {'4', '8', '11', '17', '19'}};
% wheel_number = 3;
% wheel_number = 4;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
flag = exist('folder_travel', 'var') & ...
       exist('folder_area', 'var') & ...
       exist('wheel_number', 'var');

if ~flag
    disp(['Ошибка! Проверьте наличие в рабочем пространстве ' ...
        'переменных: folder_travel, folder_area, wheel_number.']);
    return
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Основная логика работы программы

% Путь к основному скрипту и данным для оценки передаточных функций
main_path = fileparts(mfilename('fullpath'));
data_tf_path = fullfile(main_path, '..', 'data', 'tf evaluation data');
tf_path = fullfile(main_path, '..', 'data', 'tf system');

% Список файлов в директории где хранятся уже полученные данные
ceil_name_data = setdiff({dir(data_tf_path).name}, {'.', '..'});
ceil_name_tf = setdiff({dir(tf_path).name}, {'.', '..'});

fprintf('Запуск программы...\n')
disp("Список сохраненных файлов данных: " + strjoin(ceil_name_data, ', '))
fprintf(['Если в списке нет нужных файлов, либо были внесенны ...' ...
         'изменения в folder_travel или folder_area, следует ' ...
         'перезаписать данные.\n'])

prompt = 'Выбор варианта действий.';
option1 = 'Загрузить входные и выходные данные с постоянной памяти.';
option2 = 'Укомплектовать входные и выходные данные системы и сохранить.';
option3 = 'Прервать выполнения программы.';

choice = menu(prompt, option1, option2, option3);

switch choice
    case 1
        choice_user = true;
    case 2
        choice_user = false;
    case 3
        fprintf('Запуск прерван.\n')
        return
end

% Проверка наличия сохраненных данных на жестком диске
namefiles = {sprintf('data_input_outputs_signals%d.mat', wheel_number), ...
             sprintf('Ts%d.mat', wheel_number)};

flag = all(ismember(namefiles, ceil_name_data));

if flag && choice_user

    load(fullfile(data_tf_path, ...
        sprintf('data_input_outputs_signals%d.mat', wheel_number)))
    load(fullfile(data_tf_path, sprintf('Ts%d.mat', wheel_number)))

    fprintf('Данные системы успешно загруженны.\n')

else 

    if ~flag

        fprintf(['Загрузка невозможна, одного или нескольких файлов ' ...
            'нет, данные будут обновлены.\n'])
    end

    % Запуск функции создания объектов данных
    [data_input_outputs_signals, Ts] ...
    = getTotalArrays(folder_travel, folder_area, wheel_number);

    save(fullfile(data_tf_path, ...
         sprintf('data_input_outputs_signals%d.mat', wheel_number)), ...
                 'data_input_outputs_signals')
    save(fullfile(data_tf_path, ...
         sprintf('Ts%d.mat', wheel_number)), 'Ts')

    fprintf('Данные системы успешно записанны в постоянную память.\n')

end

% Создание структур данных
data_axis = iddata(data_input_outputs_signals(:, 2), ...
                   data_input_outputs_signals(:, 1), Ts);
data_rail = iddata(data_input_outputs_signals(:, 3), ...
                   data_input_outputs_signals(:, 1), Ts);
data_diag = iddata(data_input_outputs_signals(:, 4), ...
                   data_input_outputs_signals(:, 1), Ts);

clear option1 option2 option3 prompt choice_user choice flag namefiles ...
    ceil_name_data folder_area folder_travel
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Оценка передаточных функций
fprintf('Запуск оценки передаточных функций...\n')
disp("Список сохраненных файлов данных: " + strjoin(ceil_name_tf, ', '))
fprintf(['Если в списке нет нужных файлов, либо были внесенны ' ...
    'изменения в folder_travel или folder_area, следует оценить ' ...
    'передаточные \nфункция повторно. Примечание, это длительный ' ...
    'и ресурсозатратный процесс, так что если изменений не было, ' ...
    'то следует \nзагрузить данные из постоянной памяти.\n'])

prompt = 'Выбор варианта действий.';
option1 = 'Загрузить передаточные функции с постоянной памяти.';
option2 = 'Оценить передаточные функции и сохранить их.';
option3 = 'Прервать выполнения программы.';

choice = menu(prompt, option1, option2, option3);

switch choice
    case 1
        choice_user = true;
    case 2
        choice_user = false;
    case 3
        return
end

% Проверка наличия сохраненных данных на жестком диске
namefiles = {sprintf('sysd_axis%d.mat', wheel_number), ...
             sprintf('sysd_rail%d.mat', wheel_number), ...
             sprintf('sysd_diag%d.mat', wheel_number)};
    
flag = all(ismember(namefiles, ceil_name_tf));

if flag && choice_user

    load(fullfile(tf_path, sprintf('sysd_axis%d.mat', wheel_number)))
    load(fullfile(tf_path, sprintf('sysd_rail%d.mat', wheel_number)))
    load(fullfile(tf_path, sprintf('sysd_diag%d.mat', wheel_number)))

    fprintf('Данные системы успешно загруженны.\n')

else 

    if ~flag

        fprintf(['Загрузка невозможна, одного или нескольких файлов ' ...
            'нет, данные будут обновлены.\n'])

    end
        
    % Оценка передаточной функции axis
    fprintf(['Запуск оценки передаточной функции axis...\n' ...
            'Пройденные порядки модели: '])
    [sysd_axis] = evalSysTF(data_axis);
    save(fullfile(tf_path, sprintf('sysd_axis%d.mat', wheel_number)), ...
        "sysd_axis")

    % Оценка передаточной функции axis
    fprintf(['Запуск оценки передаточной функции rail...\n' ...
            'Пройденные порядки модели: '])
    [sysd_rail] = evalSysTF(data_rail);
    save(fullfile(tf_path, sprintf('sysd_rail%d.mat', wheel_number)), ...
        "sysd_rail")

    % Оценка передаточной функции axis
    fprintf(['Запуск оценки передаточной функции diag...\n' ...
            'Пройденные порядки модели: '])
    [sysd_diag] = evalSysTF(data_diag);
    save(fullfile(tf_path, sprintf('sysd_diag%d.mat', wheel_number)), ...
        "sysd_diag")

    fprintf('Данные системы успешно записанны в постоянную память.\n')

end

fprintf('Время работы составило %d минут.\n', ceil(toc / 60))

clear ceil_name_tf choice choice_user flag option1 option2 option3 ...
    prompt tf_path namefiles data_tf_path
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Сравнение данных полученных с помощью дискретной модели и реальных данных
%%
% Графическое сравнение дискретной модели sysd_axis и реальных данных
figure('Position', [100, 100, 700, 500])
% figure('Position', [100, 100, 350, 250])
compare(data_axis, sysd_axis);
ylim([-3 3])
ylabel('Нормированные ускорения, м^{-1}')
% xlabel('Время, с')
xlabel('')
% grid on
% title("axis")
title("")

% Графическое сравнение дискретной модели sysd_axis и реальных данных
figure('Position', [100, 100, 700, 500])
% figure('Position', [100, 100, 350, 250])
compare(data_rail, sysd_rail);
ylim([-3 3])
ylabel('Нормированные ускорения, м^{-1}')
% xlabel('Время, с')
xlabel('')
% grid on
% title("rail")
title("")

% Графическое сравнение дискретной модели sysd_axis и реальных данных
figure('Position', [100, 100, 700, 500])
% figure('Position', [100, 100, 350, 250])
compare(data_diag, sysd_diag);
ylim([-3 3])
ylabel('Нормированные ускорения, м^{-1}')
% xlabel('Время, с')
xlabel('')
% grid on
% title("diag")
title("")
% %%
% figure
% bode(sysd_axis)
% 
% figure
% bode(sysd_rail)
% 
% figure
% bode(sysd_diag)
% %%
% t = 0:Ts:(length(data_input_outputs_signals) - 1) * Ts;
% axis_ist = data_axis.OutputData;
% axis_restored = lsim(sysd_axis, data_axis.InputData, t);
% 
% a = 25992;
% b = 28053;
% 
% s = 50;
% [R, ~]=xcorr(axis_ist(a:b), axis_restored(a:b), 1*s, 'normalized');
% disp(max(R))
% 
% figure(4)
% plot(axis_ist(a:b), "k")
% hold on
% plot(axis_restored(a:b), "r")
% legend("Исходный S_{axis}", "Модель sysd\_axis: " + 0.75)
% ylim([-3 3])
% xlim([0 length(axis_ist(a:b))])
% ylabel('Нормированные ускорения, м^{-1}')
% xlabel('Отсчеты')
% hold off
% 
% 
% 
% rail_ist = data_rail.OutputData;
% rail_restored = lsim(sysd_rail, data_rail.InputData, t);
% 
% [R, ~]=xcorr(rail_ist(a:b), rail_restored(a:b), 1*s, 'normalized');
% disp(max(R))
% 
% figure(5)
% plot(rail_ist(a:b), "k")
% hold on
% plot(rail_restored(a:b), "r")
% legend("Исходный S_{rail}", "Модель sysd\_rail 0.5")
% ylim([-3 3])
% xlim([0 length(rail_ist(a:b))])
% ylabel('Нормированные ускорения, м^{-1}')
% xlabel('Отсчеты')
% hold off
% 
% 
% diag_ist = data_diag.OutputData;
% diag_restored = lsim(sysd_diag, data_diag.InputData, t);
% 
% [R, ~]=xcorr(diag_ist(a:b), diag_restored(a:b), 1*s, 'normalized');
% disp(max(R))
% 
% figure(6)
% plot(diag_ist(a:b), "k")
% hold on
% plot(diag_restored(a:b), "r")
% legend("Исходный S_{diag}", "Модель sysd\_diag 0.24")
% ylim([-3 3])
% xlim([0 length(diag_ist(a:b))])
% ylabel('Нормированные ускорения, м^{-1}')
% xlabel('Отсчеты')
% hold off
