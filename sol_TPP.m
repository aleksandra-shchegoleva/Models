%% функция построения графика TPP модели по заданным параметрам 
% N - время моделирования
% x - начальные значения популяции, формат: [P(0), Z(0)]
% r, K, a, beta, mu, gamma, theta - коэффициенты
% по умолчанию все аргументы равны значениям из таблицы 5, P(0) = 20, Z(0)
% = 5, N = 400

function [x] = sol_TPP(N, x, theta, r, K, a, beta, mu, gamma)

if nargin == 0
    r = 2;
    K = 108;
    a = .7;
    beta = .6;
    mu = .012;
    gamma = 5.7;
    theta = 0.5;
    N = 400;
    x = [20, 5];
end

%% Метод ode45 (Рунге-Кутта)
[M, x] = ode45(@(t, x) func(t, x, r, K, a, beta, mu, gamma, theta), [0 N], x(1,:));

x1_str = texlabel('x_1(0) = ');
x2_str = texlabel('x_2(0) = ');
r_str = texlabel('r = ');
a_str = texlabel('a = ');
b_str = texlabel('beta = ');
mu_str = texlabel('mu = ');
theta_str = texlabel('theta = ');
K_str = texlabel('K = ');
g_str = texlabel('gamma = ');

str = {strcat(x1_str,string(x(1,1))),strcat(x2_str,string(x(1,2))),strcat(r_str,string(r)),strcat(a_str,string(a)),strcat(b_str,string(beta)),strcat(mu_str,string(mu)),strcat(theta_str,string(theta)),strcat(K_str,string(K)),strcat(g_str,string(gamma))};

width  = 10;
height = 7.5;

figure('Units','inches','Position',[1 1 width height])
hold on;
box on;
grid on;
plot(M, x(:,1), 'g','Linewidth',3);
plot(M, x(:,2), 'r', 'Linewidth',3);
text(1,1,str,'Units','normalized','HorizontalAlignment', 'right', 'VerticalAlignment', 'top','FontSize',15)
axis([0 N -inf inf]);
xlabel('Время, дни');
ylabel('Популяция, ед/л');
legend({'Фитопланктон', 'Зоопланктон'}, 'Location','southeast');
ax = gca;
ax.FontSize = 20; 

figure('Units','inches','Position',[1 1 width height])
hold on;
box on;
grid on;
plot(x(:,1), x(:,2), 'b','Linewidth',3);
xlabel('x_{1}');
ylabel('x_{2}');
ax = gca;
ax.FontSize = 20;

function DXDT = func(~, x, r, K, a, beta, mu, gamma, theta)
    %функции системы ДУ
    f1 = r*x(1)*(1 - x(1)/K) - a*f(x(1),gamma)*x(2);
    f2 = beta*f(x(1),gamma)*x(2) - mu*x(2) - theta*g(x(1),gamma)*x(2);

    DXDT = [f1; f2];
end

% трофическая функция f(P)
function res = f(x, gamma)
   res = x^2/(gamma^2 + x^2);
end

% трофическая функция g(P)
function res = g(x, gamma)
   res = x^2/(gamma^2 + x^2);
end
end