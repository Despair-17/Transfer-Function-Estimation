function energyRedistribution(folder_path, main_path, wheel_choice, ...
                                                        sys_struct)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Детальное описание функции:
%   energyRedistribution - функция производит учет перераспределния энергии
%   в системе тележки вагона поезда, функция ничего не возвращет, она сох- 
%   раняет данные с поправкой под названием accel_corrected в папку с сыры-
%   ми данными. Функция рассчитвает перераспределение энергии с помощью ра-
%   нее полученных передаточных функций системы с помощью скрипта tf_evalu-
%   ation_main
%   
% Входные данные:
% - folder_path - путь к папке c сырыми файлами показаний акселерометров, 
%   которые установлены на буксах тележки вагона.
% - main_path - путь к папке main, где находятся остальные объекты (другие
%   используемые скрипты, функции и фильтр).
% - wheel_choice - выбор пользователя, что необходимо учесть.
% - sys_struct - структура содержащая все передаточные функции.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Ввод данных для учета перераспределения энергии

    % Загрузка и обработка данных
    [accel, V1mean, V2mean, Ts] = ...
                              loadAndPreprocess(folder_path, main_path, 0);

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
    % Определение индексов для учета

    % Задают длинную исследуемого участка
    len_before = 500; 
    len_after = 1499;

    idx_log = zeros(4, 2);
    idx_log(:, 2) = wheel_choice;
    [~, idx_log(1, 1)] = max(accel_filt(:, 1));
    [~, idx_log(2, 1)] = max(accel_filt(:, 2));
    [~, idx_log(3, 1)] = max(accel_filt(:, 3));
    [~, idx_log(4, 1)] = max(accel_filt(:, 4));
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % Учет перереспределения энергии
    t = 0:Ts:(len_before + len_after) * Ts;
    accel_restored = zeros(size(accel_filt));
    right_rail = [3 2 4; 4 1 3; 1 4 2; 2 3 1]';

    % Цикл учета перераспределяемой энергии
    for i = 1:4

        % Проверка необходимости учета
        log = idx_log(i, 2);
        if ~log
            continue
        end

        ind = idx_log(i, 1);
        
        % Индексовый диапазон учета
        window = ind-len_before:ind + len_after;

        % Расчет осевого перераспределения
        accel_restored(window, right_rail(1, i)) = ...
            accel_restored(window, right_rail(1, i)) + ...
            lsim(sys_struct.(sprintf('sysd_axis%d', i)), ...
                 accel_filt(window, i), t);

        % Расчет рельсового перераспределения 
        accel_restored(window, right_rail(2, i)) = ...
            accel_restored(window, right_rail(1, i)) + ...
            lsim(sys_struct.(sprintf('sysd_rail%d', i)), ...
                 accel_filt(window, i), t);
        
        % Расчет диагонального передаспределения
        accel_restored(window, right_rail(3, i)) = ...
            accel_restored(window, right_rail(1, i)) + ...
            lsim(sys_struct.(sprintf('sysd_diag%d', i)), ...
                 accel_filt(window, i), t);
        
    end

    % Учет перераспределения энергии
    accel_corrected = accel - accel_restored;

    % Возвращение от нормированных ускорений к исходным
    accel_corrected(:, 1:2) = accel_corrected(:, 1:2) * (V1mean^2);
    accel_corrected(:, 3:4) = accel_corrected(:, 3:4) * (V2mean^2);
    
    save(fullfile(folder_path, 'accel_corrected.mat'), "accel_corrected")

% Код функции заканчивается здесь
end





