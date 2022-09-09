
function [M,u,x] = sol_add_NPZ_rho_d(rho, q, T1, h, plotting, N, x, coeffs)
%     Выполняет численное моделирование NPZ-модели с аддитивным управлением
%     без расширения фазового пространства
%     Цель: достижение популяциями фитопланктона и зоопланктона
%     пропорциональных значений
%     
%     X = SOL_ADD_NPZ_RHO_D(RHO, Q, PLOTTING, N, X, COEFFS)
%     
%     ARGUMENTS
%     rho [float] коэффициент пропорциональности
%     q [float] суммарная емкость среды
%     T1 [float] параметр T1, влияющий на скорость переходного процесса
%     h [float] шаг сетки для численного метода Эйлера
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

    if nargin < 6
        x = [0.35, 0.4, .25];
        coeffs = [.05, 1, 25.003, 1.8];
    end
%     сетка для метода Эйлера
    M = 0:h:N;
    u = 0;

%     решение системы ОДУ
    [a, b, c, d] = deal(coeffs(1), coeffs(2), coeffs(3), coeffs(4));
    if q > b/c && q > (b*c + a*d*rho)/(c*d)
        x1c = (d*q - b)/(d*rho);
        x2c = b/d;
        x3c = (c*d*q - b*c - a*d*rho)/(d^2 * rho);
        fprintf('E0 q > b*c and q > (b*c + a*d*rho)/(c*d)\n');
        fprintf('%6.2f > %6.2f and %6.2f > %6.2f\n', q, b/c, q, (b*c + a*d*rho)/(c*d));
    elseif q > a*rho/c
        x1c = a/c;
        x2c = (c*q - a*rho)/c;
        x3c = 0;
        fprintf('E1 q > a*rho/c\n');
        fprintf('%6.2f > %6.2f\n', q, a*rho/c);
    else
        x1c = q/rho;
        x2c = 0;
        x3c = 0;
        fprintf('E2\n');
    end

    for i=1:length(M)-1
        f2 = c*x(i,1)*x(i,2)-d*x(i,2)*x(i,3)-a*x(i,2);
        f3 = d*x(i,2)*x(i,3) - b*x(i,3);
        psi = x(i,2) + rho*x(i,1) - q;
        u(i+1) = -psi/(T1*rho) - f2/rho - a*x(i,2) - b*x(i,3) + c*x(i,1)*x(i,2);
        f1 = a*x(i,2) + b*x(i,3) - c*x(i,1)*x(i,2) + u(i);
        
        x(i+1,1) = x(i,1) + h*f1;
        x(i+1,2) = x(i,2) + h*f2;
        x(i+1,3) = x(i,3) + h*f3;
    end
%     построение графика
    if plotting
        stpoints = [x1c, x2c, x3c];
        x1_str = texlabel('x_1(0) = ');
        x2_str = texlabel('x_2(0) = ');
        x3_str = texlabel('x_3(0) = ');
        rho_str = texlabel('rho = ');
        q_str = texlabel('q = ');
        a_str = texlabel('a = ');
        b_str = texlabel('b = ');
        c_str = texlabel('c = ');
        d_str = texlabel('d = ');
        T1_str = texlabel('T_{1} = ');
        h_str = texlabel('h = ');

        str = {strcat(x1_str,string(x(1,1))),strcat(x2_str,string(x(1,2))),...
            strcat(x3_str,string(x(1,3))),strcat(rho_str,string(rho)),...
            strcat(q_str,string(q)),strcat(a_str,string(coeffs(1))),...
            strcat(b_str,string(coeffs(2))),strcat(c_str,string(coeffs(3))),...
            strcat(d_str,string(coeffs(4))),strcat(T1_str,string(T1)),...
            strcat(h_str,string(h)), 'Метод Эйлера'};
        width      = 10;
        height     = 7.5;
        fontsize = 20;

        fig = figure('Units','inches','Position',[1 1 width height],'PaperPositionMode','auto');
        box on;
        grid on;
        hold on;
        plot(M, u, 'k','Linewidth',3);
        axis([0 N -inf inf]);
        xlabel('Время, дни');
        ylabel('Управление');
        set(fig.Children, 'FontName','Times', 'FontSize',fontsize);
        set(gca,'LooseInset', max(get(gca,'TightInset'), 0.02))
        
        fig = figure('Units','inches','Position',[1 1 width height],'PaperPositionMode','auto');
        text(1,1,str,'Units','normalized','HorizontalAlignment', 'right', 'VerticalAlignment',...
            'top','FontSize',15);
        box on;
        grid on;
        hold on;
        plot(M, x(:,1), 'g','Linewidth',3);
        axis([0 N -inf inf]);
        plot(M, x(:,2), 'r', 'Linewidth',3);
        plot(M, x(:,3), 'b', 'Linewidth',3);
        plot(M, stpoints .* ones(length(M), 3), 'k--', 'Linewidth',2);
        xlabel('Время, дни');
        ylabel('Популяция, ед/л');
        legend({'Питание (x_{1})', 'Фитопланктон (x_{2})', 'Зоопланктон (x_{3})', 'x_{1s}', 'x_{2s}', 'x_{3s}'}, 'Location','southeast');
        set(fig.Children, 'FontName','Times', 'FontSize',fontsize);
        set(gca,'LooseInset', max(get(gca,'TightInset'), 0.02))
    end

end