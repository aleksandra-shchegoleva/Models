width  = 10;
height = 7.5;
fontsize = 20;

r = 2;
K = 108;
a = .7;
beta = .6;
mu = .3654;
gamma = 5;
theta = 0.2;
T = 1;
xc = (gamma*mu)/(beta - mu - theta);
fprintf('xc = %6.2f\n', xc);

fig = figure('Units','inches','Position',[1 1 width height],'PaperPositionMode','auto');
hold on;
box on;
grid on;
x1 = linspace(0,60,5);
x2 = linspace(0,60,5);
[x1, x2] = meshgrid(x1,x2);
psi = x1 - xc;
u1 = -psi./T - r.*x1.*(1 - 1./K) + a.*x1./(gamma + x1).*x2;
x1dot = r.*x1.*(1 - x1./K) - a.*(x1./(gamma + x1)).*x2 + u1;
x2dot = beta.*(x1./(gamma + x1)).*x2 - mu.*x2 - theta.*(x1./(gamma + x1)).*x2;
streamslice(x1,x2,x1dot,x2dot);
xlabel('x_{1}');
ylabel('x_{2}');
set(fig.Children, 'FontName','Times', 'FontSize',fontsize);
set(gca,'LooseInset', max(get(gca,'TightInset'), 0.02))