function x = sol_NPZ_base(plotting, N, x, coeffs)
%     Выполняет численное моделирование базовой NPZ-модели
%     
%     X = SOL_NPZ_BASE(PLOTTING, N, X, COEFFS)
%     
%     ARGUMENTS
%     plotting [boolean] параметр для построения графика. 
%     Если принимает значение true, то график выводится в figure, 
%     если false, то возвращается только численное решение.
%     N [int] время для построения модели
%     x [1, 3] начальные условия для дифференциальных уравнений [x1(0), x2(0), x3(0)]
%     coeffs [1, 4] коэффициенты модели [a, b, c, d]
%     
%     По умолчанию 
%     x = [0.35, 0.4, .25]
%     coeffs = [.05, 1, 25.003, 1.8]
%     
%     OUTPUT
%     x [N, 3] массив с численным решением дифференциального уравнения [x1, x2, x3]
%      

    if nargin < 4
        x = [0.35, 0.4, .25];
        coeffs = [.05, 1, 25.003, 1.8];
    end
%     общая концентрация питания A = x1(t) + x2(t) + x3(t)
    A = sum(x);
    
    % вычисление стационарных точек
    [a, b, c, d] = deal(coeffs(1), coeffs(2), coeffs(3), coeffs(4));
    if c < a/A && A*d > b
        stpoints = [A, 0, 0];
        name_stpoint = 'Точка E_{0}';
    elseif A*d < b
        stpoints = [0, 0, A];
        name_stpoint = 'Точка E_{1}';
    elseif (a/A < c) && (c < a*d/(A*d - b))
        stpoints = [a/c, A - a/c, 0];
        name_stpoint = 'Точка E_{2}';
    else
        stpoints = [(a - b + A*d)/(c + d), b/d, (A*c*d - a*d - b*c)/(d^2 + c*d)];
        name_stpoint = 'Точка E_{3}';
    end
    fprintf(name_stpoint);
    fprintf('\nx_1s = %6.2f, x_2s = %6.2f, x_3s = %6.2f\n',stpoints);

%     решение системы ОДУ
    [M, x] = ode45(@(t, x) odenpz(t, x), [0 N], x(1,:));
    
%     построение графика
    if plotting
        x1_str = texlabel('x_1(0) = ');
        x2_str = texlabel('x_2(0) = ');
        x3_str = texlabel('x_3(0) = ');
        a_str = texlabel('a = ');
        b_str = texlabel('b = ');
        c_str = texlabel('c = ');
        d_str = texlabel('d = ');

        str = {strcat(x1_str,string(x(1,1))),strcat(x2_str,string(x(1,2))),...
            strcat(x3_str,string(x(1,3))),strcat(a_str,string(coeffs(1))),...
            strcat(b_str,string(coeffs(2))),strcat(c_str,string(coeffs(3))),...
            strcat(d_str,string(coeffs(4)))};
        width      = 10;
        height     = 7.5;
        fontsize = 20;

        fig = figure('Units','inches','Position',[1 1 width height],'PaperPositionMode','auto');
        text(1,1,str,'Units','normalized','HorizontalAlignment', 'right', 'VerticalAlignment',...
            'top','FontSize',15);
        box on;
        grid on;
        hold on;
        plot(M, x(:,1), 'g','Linewidth',3);
        axis([0 N -inf A+0.1]);
        plot(M, x(:,2), 'r', 'Linewidth',3);
        plot(M, x(:,3), 'b', 'Linewidth',3);
        plot(M, stpoints .* ones(length(M), 3), 'k--', 'Linewidth',2);
        xlabel('Время, дни');
        ylabel('Популяция, ед/л');
        legend({'Питание (x_{1})', 'Фитопланктон (x_{2})', 'Зоопланктон (x_{3})', 'x_{1s}', 'x_{2s}', 'x_{3s}'}, 'Location','best');
        set(fig.Children, 'FontName','Times', 'FontSize',fontsize);
        set(gca,'LooseInset', max(get(gca,'TightInset'), 0.02))
    end
    function [output] = odenpz(~, x)
       f1 = a*x(2) + b*x(3) - c*x(1)*x(2);
       f2 = c*x(1)*x(2) - d*x(2)*x(3) - a*x(2);
       f3 = d*x(2)*x(3) - b*x(3);
       
       output = [f1; f2; f3];
    end
end