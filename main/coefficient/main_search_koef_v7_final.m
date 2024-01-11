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

% Сглаживание сигнала по 8 точкам
% n = 8;
% accel_1 = smoothdata(accel_1, 'movmean', n);
% accel_2 = smoothdata(accel_2, 'movmean', n);
% accel_3 = smoothdata(accel_3, 'movmean', n);
% accel_4 = smoothdata(accel_4, 'movmean', n);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Общий графики участка ускорений акселерометров

figure(1)
plot(1:length(accel_1), accel_1,'k', 1:length(accel_2), accel_2, 'r', ...
    1:length(accel_3), accel_3, 'b', 1:length(accel_4), accel_4,'g')
    xlabel('Отсчеты')
    ylabel('Нормированные ускорения, м^{-1}')
%     ylabel('Ускорения, м/с^2')
%     legend('1 Следом идущее','2 Впереди идущее', ...
%         '3 Следом идущее','4 Впереди идущее')
    legend('1 Впереди идущее','2 Следом идущее', ...
           '3 Впереди идущее','4 Следом идущее')
    legend('1 Leading','2 Trailing', '3 Leading', '4 Trailing')
    grid on
    xlabel('sample #')
    ylabel('adjusted acceleration, m⁻¹')
%     set(gca, 'XTick', 0:250:length(accel_1))
%     xlim([2250 4250])
%     ylim([-6 12])
    xlim([0 length(accel_1)])
    ylim([min(min(accel_1) - 4, min(accel_3) - 4) ...
        max(max(accel_1) + 4, max(accel_3) + 4)])
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
n = 500; k = 1499;

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

% Графики совмещенных ускорений
figure(2)
plot(1:length(accel_1_comb), accel_1_comb,'k', ...
    1:length(accel_2_comb), accel_2_comb, 'r')
    xlabel('Отсчеты')
    ylabel('Нормированные ускорения, м^{-1}')
%     ylabel('Ускорения, м/с^2')
    legend('1 Следом идущее','2 Впереди идущее')
    grid on  
    set(gca, 'XTick', 0:250:length(accel_1_comb))
    ylim([min(min(accel_1_comb) - 1, min(accel_3_comb) - 1) ...
        max(max(accel_1_comb) + 1, max(accel_3_comb) + 1)])

figure(3)
plot(1:length(accel_3_comb), accel_3_comb, 'b', ...
    1:length(accel_4_comb), accel_4_comb,'g')
    xlabel('Отсчеты')
    ylabel('Нормированные ускорения, м^{-1}')
%     ylabel('Ускорения, м/с^2')
    legend('3 Следом идущее','4 Впереди идущее')
    grid on
    set(gca, 'XTick', 0:250:length(accel_3_comb))
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
%% Рассчет коэффициентов 

% Быстрое преобразование Фурье                  
L = length(accel_1_comb);
Fs = 31500;

Res_K = zeros(4,5);
for i = 1:1:4
    
    if i == 1
        m = index_1;
    elseif i == 2
        m = index_2;
    elseif i == 3
        m = index_3;
    elseif i == 4
        m = index_4;
    end  

    % Преобразование Фурье для модуля 1 
    Y1 = fft(accel_1(m - n:m + k));                   
    P2_1 = abs(Y1 / L);
    P1_1 = P2_1(1:round(L / 2) + 1);
    P1_1(2:end-1) = 2 * P1_1(2:end - 1);

    % Преобразование Фурье для модуля 2
    Y2 = fft(accel_2(m - n:m + k));                   
    P2_2 = abs(Y2 / L);
    P1_2 = P2_2(1:round(L / 2) + 1);
    P1_2(2:end-1) = 2 * P1_2(2:end - 1);

    % Преобразование Фурье для модуля 3
    Y3 = fft(accel_3(m - n:m + k));                   
    P2_3 = abs(Y3 / L);
    P1_3 = P2_3(1:round(L / 2) + 1);
    P1_3(2:end-1) = 2 * P1_3(2:end - 1);

    % Преобразование Фурье для модуля 4
    Y4 = fft(accel_4(m - n:m + k));                   
    P2_4 = abs(Y4 / L);
    P1_4 = P2_4(1:round(L / 2) + 1);
    P1_4(2:end-1) = 2 * P1_4(2:end - 1);

    clear Y1 Y2 Y3 Y4 P2_1 P2_2 P2_3 P2_4 T t

    % Частотный диапазон
    f = Fs*(0:(L/2)) / L;

    % Графики амлитудного спектров
    % b = 500;
    b = length(f);
    
    o = 2;
    P1_1 = smoothdata(P1_1, 'sgolay', o);
    P1_2 = smoothdata(P1_2, 'sgolay', o);
    P1_3 = smoothdata(P1_3, 'sgolay', o);
    P1_4 = smoothdata(P1_4, 'sgolay', o);

%     figure('Position', [100, 100, 700, 500])
    figure('Position', [100, 100, 350, 250])
    plot(f(1:b), P1_1(1:b), 'k', f(1:b), P1_2(1:b), 'r', f(1:b), ...
        P1_3(1:b), 'b', f(1:b), P1_4(1:b), 'g') 
%         title('Амплитудный спектр')
        xlabel('f, Гц')
%         ylabel('Амплитуда, м^{-1}')
        ylabel('amplitude, м^{-1}')
        grid on
%         legend('1 Следом идущее', '2 Впереди идущее', ...
%                '3 Следом идущее', '4 Впереди идущее')
        legend('1 Leading','2 Trailing', '3 Leading', '4 Trailing')
%         xlim([0 500])

    clear f
    % Максимальная амплитуда спектров
    b = 300;
    Amax_1 = max(P1_1(1:b));
    Amax_2 = max(P1_2(1:b));
    Amax_3 = max(P1_3(1:b));
    Amax_4 = max(P1_4(1:b));

    % От колеса 1 к остальным
    if max([Amax_1, Amax_2, Amax_3, Amax_4]) == Amax_1
        K_axis = roundn(Amax_3 / Amax_1, -3);
        K_rail = roundn(Amax_2 / Amax_1, -3);
        K_diag = roundn(Amax_4 / Amax_1, -3);
        Res_K(i,1) = K_axis; Res_K(i,2) = K_rail; Res_K(i,3) = K_diag; 
        Res_K(i,4) = Amax_1; Res_K(i,5) = str2double(lot_number);
    
    % От колеса 2 к остальным    
    elseif max([Amax_1, Amax_2, Amax_3, Amax_4]) == Amax_2
        K_axis = roundn(Amax_4 / Amax_2, -3);
        K_rail = roundn(Amax_1 / Amax_2, -3);
        K_diag = roundn(Amax_3 / Amax_2, -3);
        Res_K(i,1) = K_axis; Res_K(i,2) = K_rail; Res_K(i,3) = K_diag; 
        Res_K(i,4) = Amax_2; Res_K(i,5) = str2double(lot_number);
    
    % От колеса 3 к остальным  
    elseif max([Amax_1, Amax_2, Amax_3, Amax_4]) == Amax_3
        K_axis = roundn(Amax_1 / Amax_3, -3);
        K_rail = roundn(Amax_4 / Amax_3, -3);
        K_diag = roundn(Amax_2 / Amax_3, -3);   
        Res_K(i,1) = K_axis; Res_K(i,2) = K_rail; Res_K(i,3) = K_diag; 
        Res_K(i,4) = Amax_3; Res_K(i,5) = str2double(lot_number);
    
    % От колеса 4 к остальным    
    elseif max([Amax_1, Amax_2, Amax_3, Amax_4]) == Amax_4
        K_axis = roundn(Amax_2 / Amax_4, -3);
        K_rail = roundn(Amax_3 / Amax_4, -3);
        K_diag = roundn(Amax_1 / Amax_4, -3); 
        Res_K(i,1) = K_axis; Res_K(i,2) = K_rail; Res_K(i,3) = K_diag; 
        Res_K(i,4) = Amax_4; Res_K(i,5) = str2double(lot_number);
    
    end
end
clear P1_1 P1_2 P1_3 P1_4 m b i Amax_1 Amax_2 Amax_3 Amax_4 L n k ...
    K_axis K_diag K_rail koef

fprintf("Амплитудные коэффициенты при проезде %s на участке %s.\n", ...
                                                  road_name, lot_number)
disp(Res_K)

% close(4:8)
%%
% Выбор длины участка
n = 500; k = 1499;

figure
% plot(1:length(accel_1), accel_1,'k', 1:length(accel_2), accel_2, 'r', ...
%     1:length(accel_3), accel_3, 'b', 1:length(accel_4), accel_4,'g')

ind = 2749;
% 
% signal_koef = accel_2(ind - n:ind + k) * 0.28;
% signal = accel_4(ind - n:ind + k);
% 
signal_koef = accel_1(ind - n:ind + k) * 0.28;
signal = accel_3(ind - n:ind + k);

plot(1:length(signal_koef), signal_koef, 'k', 1:length(signal), signal,'r')
    xlabel('Отсчеты')
    xlabel('sample #')
    ylabel('Нормированные ускорения, м^{-1}')
    ylabel('adjusted acceleration, m⁻¹')
    legend('Корректировка', 'Исходный')
    legend('Adjustment', 'Raw')
%     ylabel('Ускорения, м/с^2')
%     legend('1 Следом идущее','2 Впереди идущее', ...
%         '3 Следом идущее','4 Впереди идущее')
%     legend('1 Впереди идущее','2 Следом идущее', ...
%            '3 Впереди идущее','4 Следом идущее')
    grid on
%     set(gca, 'XTick', 0:250:length(accel_1))
%     xlim([2250 4250])
%     ylim([-4 4])
%     xlim([0 length(accel_1)])
%     ylim([min(min(accel_1) - 4, min(accel_3) - 4) ...
%         max(max(accel_1) + 4, max(accel_3) + 4)])



