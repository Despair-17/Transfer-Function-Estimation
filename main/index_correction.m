function [correction_2, correction_4] = ...
          index_correction(accel_1, accel_2, accel_3, accel_4)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Детальное описание функции:
%   index_correction - функция рассчитывает поправки к грубо оцененным ин-
%   дексам с помощью взаимной корреляции.
%
% Входные аргументы:
%   accel_i - сигнал i акселерометра.
%
% Выходные аргументы:
%   correction_2 - поправка индекса к второму акселерометру
%   correction_4 - поправка индекса к четвертому акселерометру
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Код функции начинается здесь
    i = 0;
    s = length(accel_1);
    b = s * 0.75;
    R1 = zeros(s,1);
    R2 = zeros(s,1);

    while i<s
        y1 = circshift(accel_1, -round(b) + i);
        y2 = accel_2(1:end);
        Ri1 = corr(y1,y2);
        R1(i+1,:) = Ri1;
    
        y3 = circshift(accel_3, -round(b) + i);
        y4 = accel_4(1:end);
        Ri2 = corr(y3,y4);
        R2(i+1,:) = Ri2;
    
        clear Ri1 Ri2 y1 y2 y3 y4
        i = i + 1;
    end

    [~,k1] = max(R1);
    correction_2 = (k1) - round(b);
    [~,k2] = max(R2);
    correction_4 = (k2) - round(b);

% Код функции заканчивается здесь
end