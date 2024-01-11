%% 
clear; close all; clc 
%%
load('koef.txt')
roef_sort = sortrows(koef, 4);

K_axis = roef_sort(:,1);
K_rail = roef_sort(:,2);
K_diag = roef_sort(:,3);

K_axis_m = mean(K_axis);
K_rail_m = mean(K_rail);
K_diag_m = mean(K_diag);

%%
% K_axis = rmoutliers(K_axis, 'mean');
% K_rail = rmoutliers(K_rail, 'mean');
% K_diag = rmoutliers(K_diag, 'mean');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure(1)
plot(roef_sort(:,4), K_axis, "k")
xlabel("Амплитуда, м^{-1}")
ylabel("Коэффициент K_{axis}")
xlim([0.2 2.2])
ylim([0 0.5])
% legend('axis')
grid on
set(gca, 'YTick', 0:0.05:1)
% set(gca, 'XTick', 0:0.5:7)

% 2
figure(2)
plot(roef_sort(:,4), K_rail, "k")
xlabel("Амплитуда, м^{-1}")
ylabel("Коэффициент K_{rail}")
xlim([0.2 2.2])
ylim([0 0.5])
% legend('rail')
grid on
set(gca, 'YTick', 0:0.05:1)
% set(gca, 'XTick', 0:0.5:7)

figure(3)
plot(roef_sort(:,4), K_diag, "k")
xlabel("Амплитуда, м^{-1}")
ylabel("Коэффициент K_{diag}")
xlim([0.2 2.2])
ylim([0 0.5])
% legend('diag')
grid on
set(gca, 'YTick', 0:0.05:1)
% set(gca, 'XTick', 0:0.5:7)


x = roef_sort(:,4);
x_a = x(1):1:x(end);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 1
p1 = polyfit(x,K_axis, 1);
y1 = polyval(p1,x_a);
figure(4)
hold on
plot(x_a, y1, "r")
xlabel("Амплитуда, м^{-1}")
ylabel("Коэффициент K_{axis}")
xlim([0.2 2.2])
ylim([0 0.6])
% plot(roef_sort(:,4), K_axis, "k")
% legend('axis')
grid on
% set(gca, 'YTick', 0:0.05:1)

% 2
p2 = polyfit(x,K_rail, 1);
y2 = polyval(p2,x_a);
figure(5)
hold on
plot(x_a, y2, "r")
xlabel("Амплитуда, м^{-1}")
ylabel("Коэффициент K_{rail}")
xlim([0.2 2.2])
ylim([0 0.6])
% plot(roef_sort(:,4), K_rail, "k")
% legend('rail')
grid on
% set(gca, 'YTick', 0:0.05:1)
% 
% 
p3 = polyfit(x, K_diag, 1);
y3 = polyval(p3,x_a);
figure(6)
hold on
plot(x_a, y3, "r")
xlabel("Амплитуда, м^{-1}")
ylabel("Коэффициент K_{diag}")
xlim([0.2 2.2])
ylim([0 0.6])
% plot(roef_sort(:,4), [K_diag; K_diag(1)], "k")
% legend('diag')
grid on
% set(gca, 'YTick', 0:0.05:1)

% K_axis_m = mean(K_axis);
% K_rail_m = mean(K_rail);
% K_diag_m = mean(K_diag);
close(4:6)


figure
line([0.3597 2.0854], [0.255202 0.268489], 'Color', 'r'); 
xlabel("Амплитуда, м^{-1}")
ylabel("Коэффициент K_{axis}")
xlim([0.2 2.2])
ylim([0 0.5])
grid on
set(gca, 'YTick', 0:0.05:1)

figure
line([0.3597 2.0854], [0.177086 0.217052], 'Color', 'r'); 
xlabel("Амплитуда, м^{-1}")
ylabel("Коэффициент K_{rail}")
xlim([0.2 2.2])
ylim([0 0.5])
grid on
set(gca, 'YTick', 0:0.05:1)

figure
line([0.3597 2.0854], [0.155968 0.207141], 'Color', 'r'); 
xlabel("Амплитуда, м^{-1}")
ylabel("Коэффициент K_{diag}")
xlim([0.2 2.2])
ylim([0 0.5])
grid on
set(gca, 'YTick', 0:0.05:1)