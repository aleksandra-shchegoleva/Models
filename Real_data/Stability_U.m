N = 100;
h = .01;
M = 0:h:N;

e1 = 10^6; %критерий стабильности
flag = 0;

alpha2 = 0.12;
beta1 = 0.13;
beta2 = 0.05;
T1 = 1;
T2 = 1;

mas = [];
flag_find = 1;
x = [4 20 0.04];
% for X1 = 1:1:30
% for X2 = 1:1:30

xc = alpha2/beta2;
e2 = .1 * xc; %критерий достижения цели
    
for a1 = 0.01:0.01:100
    x = [4 20 a1];
% for alpha2 = 0.1:0.1:1
% for beta1 = 0.01:0.01:1
% for beta2 = 0.01:0.01:1
% for T1 = 0:.1:20
% for T2 = 0:100
flag_xc = 0;
    flag = 0;
    x = solODExc(h, N, x, xc, alpha2, beta1, beta2, T1, T2, 2, false);
    
    if sum(find(x(:,1) < 0)) > 0 || isnan(x(end,1))
       flag = 1;
    end
    if sum(find(x(:,2) < 0)) > 0 || isnan(x(end,2))
       flag = 1;
    end
    if isnan(x(end,3))
       flag = 1;
    end
    if flag == 0
        for i=1:length(M)
            if abs(x(i,1) - xc) <= e2 && abs(x(i + 1,1) - xc) <= e2 && abs(x(end,1) - xc) <= e2
                mas(end + 1,:) = [x(1,1) x(1,2) T1 T2 M(i) x(1,3) alpha2 beta1 beta2];
%                 flag_xc = 1;
%                 dx2=[0 diff(x2)];
%                 [x2_1,t1]=find(dx2<0);
%                 [x2_2,t2]=find(dx2>0);
%                 [x2_3,t3]=find(dx2==0);
%                 S1 = sum(ismember(M(N/2+1:end), t1));
%                 S2 = sum(ismember(M(N/2+1:end), t2));
%                 S3 = sum(ismember(M(N/2+1:end), t3));
%                 S1_1 = sum(ismember(M(1:N/2), t2));
%                 S2_1 = sum(ismember(M(1:N/2), t1));
%                 if S1 > 0.8*(N - (N/2+1))
%                     dataLess0(end+1,:) = [x1(1) x2(1) T1 T2 i tx2 alpha2 beta1 beta2];
%                 end
%                 if S2 > 0.8*(N - (N/2+1))
%                     dataMore0(end+1,:) = [x1(1) x2(1) T1 T2 i tx2 alpha2 beta1 beta2];
%                 end
%                 if S3 > 0.8*(N - (N/2+1))
%                     dataEq0(end+1,:) = [x1(1) x2(1) T1 T2 i tx2 alpha2 beta1 beta2];
%                 end
%                 if S1_1 > 0.05*(N/2)
%                     dataMore0X2(end+1,:) = [x1(1) x2(1) T1 T2 i tx2 alpha2 beta1 beta2];
%                 end
%                 if S2_1 > 0.05*(N/2)
%                     dataLess0X2(end+1,:) = [x1(1) x2(1) T1 T2 i tx2 alpha2 beta1 beta2];
%                 end
                break;
            end
        end
    end
end
% end
% end

%% построение графика

i = 2312;
row = num2cell(mas(i,:));
x = [];
[x(1,1), x(1,2), T1, T2, t, x(1,3), alpha2, beta1, beta2] = deal(row{:});
x = solODExc(0.01, 500, [4 20 24.13], xc, alpha2, beta1, beta2, T1, T2, 2, true);

%%
plot(mas(:,6),mas(:,5),'o','LineWidth',3);
xlabel("\alpha_{1}(0)");
ylabel("t_{достижения цели}, дни");
ax = gca;
ax.FontSize = 30;

%%
clear;
t = 0:50;
a2 = 0.05;
b2 = 0.3;
d = 2;
r = 2;
x2 = (-a2/b2 - d)./(r*(exp((a2 + b2*d).*t) - 1));
plot(t,x2, 'LineWidth',3);
axis([0 inf -inf 0]);
xlabel('t');
ylabel('x_{2}(t)');
ax = gca;
ax.FontSize = 20;

%%
clear;
t = 0:50;
a2 = 0.05;
b2 = 0.3;
d = 2;
r = 2;
x2 = ((-a2/b2 + d)*exp(a2 + b2*d.*t))./(r*(exp(a2 + b2*d.*t) - exp(a2.*t + b2*d)));
plot(t,x2, 'LineWidth',3);
axis([0 inf 0 inf]);
xlabel('t');
ylabel('x_{2}(t)');
ax = gca;
ax.FontSize = 20;