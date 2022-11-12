% syms x1 x2 x3 a b c d A
% f1 = a*x2 + b*x3 - c*x1*x2 == 0;
% f2 = c*x1*x2 - d*x2*x3 - a*x2 == 0;
% f3 = d*x2*x3 - b*x3 == 0;
% f4 = x1 + x2 + x3 == A;
% 
% x = [f1, f2, f3, f4];
% x = solve(x,[x1 x2 x3]);

%% npz model additive xc
syms x1 x2 x3 a b c d T1 T2 xc
psi = x2 - xc;
phi = (-psi/T2 + d*x2*x3 + a*x2) / (c*x2);
psi1 = x1 - phi;
f2 = c*x1*x2 - d*x2*x3 - a*x2;
f3 = d*x2*x3 - b*x3;
dphidt = ((-f2/T2 + d*(f2*x3 + x2*f3) + a*f2)*x2 - (-psi/T2 + d*x2*x3 + ...
    a*x2)*f2) / (c * x2^2);
u = -psi1/T1 - a*x2 - b*x3 + c*x1*x2 + dphidt;
f1 = a*x2 + b*x3 - c*x1*x2 + u;

x = [f1, f2];
x = solve(x,[x1 x3]);

%% npz model additive rho q
syms x1 x2 x3 a b c d T1 rho q
f2 = c*x1*x2-d*x2*x3-a*x2;
f3 = d*x2*x3 - b*x3;
psi = x2 + rho*x1 - q;
u = -psi/(T1*rho) - f2/rho - a*x2 - b*x3 + c*x1*x2;
f1 = a*x2 + b*x3 - c*x1*x2 + u;

x = [f1, f2, f3];
x = solve(x,[x1 x2 x3]);

%% tpp model additive xc
syms x1 x2 r K a beta mu theta gamma xc T
psi = x1 - xc;
u = -psi/T - r*x1*(1 - x1/K) + a*(x1/(gamma + x1))*x2;
f1 = r*x1*(1 - x1/K) - a*(x1/(gamma + x1))*x2 + u;
f2 = beta*(x1/(gamma + x1)) - mu*x2 - theta*(x1/(gamma + x1))*x2;

x = [psi, f1, f2];
x = solve(x,[x1 x2]);

%% tpp model additive rho d
syms x1 x2 r K a beta mu theta gamma rho d T
psi = x1 - rho*x2 + d;
f2 = beta*(x1/(gamma + x1)) - mu*x2 - theta*(x1/(gamma + x1))*x2;
u = -psi/T + rho*f2 - r*x1*(1 - x1/K) + a*(x1/(gamma + x1))*x2;
f1 = r*x1*(1 - x1/K) - a*(x1/(gamma + x1))*x2 + u;

x = [psi, f1, f2];
x = solve(x,[x1 x2]);

%% npz big fase
syms x1 x2 x3 a b c d T1 T2 rho q
f2 = c*x1*x2-d*x2*x3-a*x2;
f3 = d*x2*x3 - b*x3;
psi = x2 + rho*x3 - q;
phi = (-psi/T2 + d*x2*x3 + a*x2 - rho*(d*x2*x3 - b*x3))/(c*x2);
dphidt = ((-(f2+rho*f3)/T2 + d*(f2*x3 + x2*f3) + a*f2 - ...
    rho*(d*(f2*x3 + x2*f3) - b*f3))*x2 - ...
    (-psi/T2 + d*x2*x3 + a*x2 - rho*(d*x2*x3 - ...
    b*x3))*f2)/(c*x2^2);
psi1 = x1 - phi;
u = -psi1/T1 - a*x2 - b*x3 + c*x1*x2 + dphidt;
f1 = a*x2 + b*x3 - c*x1*x2 + u;

x = [psi, psi1, f1, f2, f3];
x = solve(x,[x1 x2 x3]);

%% npz new model
syms x1 x2 x3 lambda1 lambda2 lambda3 a11 a22 r K mu gamma r1 r2 beta1 beta2 K1 K2

% f1 = x1*(r*(1 - x1/K) - a11*x1 - (x1/(lambda1+x1))*x2 - (x1/(lambda2+x1))*x3);
% f2 = x2*(r*(1 - x2/K) - (x2/(lambda1+x2))*x1 - a22*x2 - (x2/(lambda2+x2))*x3);
% f3 = -x3*(mu + gamma*x3 + (x2/(lambda3+x2))*x1 - (x1/(lambda2+x1))*x2);
% x = [f1, f2, f3];
% x = solve(x,[x1 x2 x3]);

% f1 = x1*(r*(1 - x1/K) - (x1/(lambda1+x1))*x2 - (x1/(lambda2+x1))*x3);
% f2 = x2*(r*(1 - x2/K) - (x2/(lambda1+x2))*x1 - (x2/(lambda2+x2))*x3);
% f3 = -x3*(mu + (x2/(lambda3+x2))*x1 - (x1/(lambda2+x1))*x2);
% x = [f1, f2, f3];
% x = solve(x,[x1 x2 x3]);

f1 = (r1*(1 - x1/K) - (x1/(lambda1+x1))*x3);
f2 = (r2*(1 - x2/K) - (x2/(lambda2+x2))*x3);
f3 = (mu + (x2/(lambda2+x2))*x1 - (x1/(lambda1+x1))*x2);

% x3_1 = solve(f1, x3);
% x3_2 = solve(f2, x3);
% 
% x1 = solve(x3_1 == x3_2, x1);
% f3 = (mu + (x2/(lambda3+x2))*x1 - (x1/(lambda1+x1))*x2);
% solve(f3, x2)

% f3_1 = (mu + (x2/(lambda3+x2))*x1(1) - (x1(1)/(lambda1+x1(1)))*x2);
% solve(f3_1, x2)

x = solve([f1, f2, f3],[x1, x2, x3])

for i=1:3
   fprintf('Точка %6.2f\n', i);
   fprintf('x1= %100s\n', x.x1(i,1));
%     pretty()
   fprintf('x2= %100s\n', x.x2(i,1));
   fprintf('x3= %100s\n\n', x.x3(i,1));
end

