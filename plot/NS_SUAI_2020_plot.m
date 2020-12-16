font = 20; %размер шрифта длЯ подписи осей (кроме рисунка 1)
font_fig1 = 16; %размер шрифта длЯ подписи оси y и названий графиков (рисунок 1)
font_fig1_x = 14; %размер шрифта длЯ подписис оси x (рисунок 1)
font_fig1_legend = 11; %размер шрифта длЯ легенды (рисунок 1)
width = 3; %ширина линий на рисунках
%% €‘•Ћ„Ќ›… „ЂЌЌ›… ‘Ћ„…ђ†€’‘џ ‚ ”Ђ‰‹Ђ• DATAx.MAT, ѓ„… х - ЌЋЊ…ђ ђ€‘“ЌЉЂ. „‹џ ђ€‘“ЌЉЂ 1 ‚‘… €‘•Ћ„Ќ›… „ЂЌЌ›…
%% “†… ‚›ЃђЂЌ›, „‹џ ђ€‘“ЌЉЂ 12 ”Ђ‰‹› €Њ…ћ’ ЌЂ‡‚ЂЌ€… DATA12-x, ѓ„… х - ђ€‘“ЌЋЉ Ђ €‹€ Ѓ
%% ‚ ”“ЌЉ–€€ LOAD ‚›Ѓ€ђ€ђЂ…’‘џ ”Ђ‰‹ „‹џ ЏЋ‘’ђЋ…Ќ€џ ѓђЂ”€ЉЂ, “ ЉЂ†„ЋѓЋ ѓђЂ”€ЉЂ “ЉЂ‡ЂЌЋ „‹џ ЉЂЉ€• ђ€‘“ЌЉЋ‚ ЋЌ ‘’ђЋ€’‘џ

%% ђ€‘“ЌЋЉ 1
% load('data1-1.mat');
% subplot(221);
% plot(M, x1, 'g','Linewidth',width);
% axis([0 N -inf inf]);
% hold on;
% plot(M, x2, 'r', 'Linewidth',width);
% hold on;
% plot(M, alpha1, 'b', 'Linewidth', width);
% hold on;
% plot(M, ones(length(M)).*xc, '--k','Linewidth', width-1);
% title(['T_1 = ',num2str(T1), ', T_2 = ', num2str(T2)], 'FontSize', font_fig1);
% xlabel('Время', 'FontSize', font_fig1_x);
% ylabel('Численность', 'FontSize', font_fig1);
% legend('Жертвы', 'Хищники','Питание', 'Целевое значение', 'FontSize', font_fig1_legend);
% 
% load('data1-2.mat');
% subplot(222);
% plot(M, x1, 'g','Linewidth',width);
% axis([0 N -inf inf]);
% hold on;
% plot(M, x2, 'r', 'Linewidth',width);
% hold on;
% plot(M, alpha1, 'b', 'Linewidth', width);
% hold on;
% plot(M, ones(length(M)).*xc, '--k','Linewidth',width-1);
% title(['T_1 = ',num2str(T1), ', T_2 = ', num2str(T2)], 'FontSize', font_fig1);
% xlabel('Время', 'FontSize', font_fig1_x);
% ylabel('Численность', 'FontSize', font_fig1);
% legend('Жертвы', 'Хищники','Питание', 'Целевое значение', 'FontSize', font_fig1_legend);
% 
% load('data1-1.mat');
% subplot(223);
% plot(M, x1, 'g','Linewidth',width);
% axis([0 N -inf inf]);
% hold on;
% plot(M, x2, 'r', 'Linewidth',width);
% hold on;
% plot(M, alpha1, 'b', 'Linewidth', width);
% hold on;
% plot(M, ones(length(M)).*xc, '--k','Linewidth',width-1);
% title(['T_1 = ',num2str(T1), ', T_2 = ', num2str(T2)], 'FontSize', font_fig1);
% xlabel('Время', 'FontSize', font_fig1_x);
% ylabel('Численность', 'FontSize', font_fig1);
% legend('Жертвы', 'Хищники','Питание', 'Целевое значение', 'FontSize', font_fig1_legend);
% 
% load('data1-3.mat');
% subplot(224);
% plot(M, x1, 'g','Linewidth',width);
% axis([0 N -inf inf]);
% hold on;
% plot(M, x2, 'r', 'Linewidth',width);
% hold on;
% plot(M, alpha1, 'b', 'Linewidth', width);
% hold on;
% plot(M, ones(length(M)).*xc, '--k','Linewidth',width-1);
% title(['T_1 = ',num2str(T1), ', T_2 = ', num2str(T2)], 'FontSize', font_fig1);
% xlabel('Время', 'FontSize', font_fig1_x);
% ylabel('Численность', 'FontSize', font_fig1);
% legend('Жертвы', 'Хищники','Питание', 'Целевое значение', 'FontSize', font_fig1_legend);


%% РИСУНКИ 2-3
% load('data3.mat');
% plot(M, x1, 'g','Linewidth',width);
% axis([0 N -inf inf]);
% hold on;
% plot(M, x2, 'r', 'Linewidth',width);
% hold on;
% plot(M, alpha1, 'b', 'Linewidth', width);
% hold on;
% plot(M, ones(length(M)).*xc, '--k','Linewidth',width-1);
% xlabel('Время', 'FontSize', font);
% ylabel('Численность', 'FontSize', font);
% legend('Жертвы', 'Хищники','Питание', 'Целевое значение', 'FontSize', font);


%% РИСУНКИ 4,7 (ФАЗОВЫЙ ПОРТРЕТ)
% load('data4.mat');
% plot(x1, x2, 'b','Linewidth',width);
% xlabel('X_1', 'FontSize', font);
% ylabel('X_2', 'FontSize', font);


%% РИСУНКИ 5,8 (ДОСТИЖЕНИЕ ЦЕЛИ)
% load('data8.mat');
% plot(M, PSI, 'g','Linewidth',width);
% xlabel('Время', 'FontSize', font);
% ylabel('{\psi}', 'FontSize', font);


%% РИСУНКИ 6,9,10,11 (СИСТЕМА С ОГРАНИЧЕНИЕМ)
% load('data11.mat');
% plot(M, x1, 'g','Linewidth',width);
% axis([0 N -inf inf]);
% hold on;
% plot(M, x2, 'r', 'Linewidth',width);
% hold on;
% plot(M, alpha1, 'b', 'Linewidth', width);
% hold on;
% plot(M, ones(length(M)).*xc, '--k','Linewidth',width-1);
% hold on;
% plot(M, ones(length(M)).*B, '--m','Linewidth',width-1);
% xlabel('Время', 'FontSize', font);
% ylabel('Численность', 'FontSize', font);
% % для рисунка 6 легенда расположена справа внизу
% % legend('Жертвы', 'Хищники','Питание', 'Целевое значение', 'Ограничение B', 'FontSize', font, 'Location', 'southeast');
% % для остальных легенда расположена справа сверху
% legend('Жертвы', 'Хищники','Питание', 'Целевое значение', 'Ограничение B', 'FontSize', font);


%% РИСУНОК 12 (3D графики)
% load('data12-1.mat');
% plot3(mas(:,1), mas(:,2), mas(:,3), 'o');
% xlabel('h', 'FontSize', font);
% ylabel('T2', 'FontSize', font);
% zlabel('T1', 'FontSize', font);