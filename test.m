clear; close all; clc;
addpath(genpath('./Lib/Linear/'))
addpath(genpath('./Lib/NonLinear/'))
A = [4 -1 0; -1 4 -1; 0 -1 4];
b = [2 4 10]';
x0 = zeros(size(b));
w = 2*(1 + sqrt(7/8))^-1;

x1 = JacobiMethod(A,b,x0);          % ./Lib/Linear/
x2 = Gauss_SeidelMethod(A,b,x0);    % ./Lib/Linear/
x3 = SORMethod(A,b,x0,w);           % ./Lib/Linear/
x4 = ConjGradMethod(A,b,x0);        % ./Lib/Linear/

f = @(x) exp(-x)-x; fp = @(x) -exp(-x)-1;
x5 = NewtonsMethod(f,fp,5.25,100); % ./Lib/NonLinear/

f = @(x) exp(x);
dt = linspace(.001,1,1000); x = 1;
dydt = FD_012_1(f,x,dt);
figure(1); hold on;
plot(log(dt),log(dydt-exp(1)));
plot(log(dt),log(dt.^2)+2 );
xlabel("$\ln(h)$",'Interpreter','latex');
ylabel("$\ln(f'-\tilde{f'})$",'Interpreter','latex');
title('Error Plot');


function dydt = FD_012_1(f,x,dt)
    dydt = (-3*f(x) + 4*f(x+dt) -f(x+2*dt))/2./dt;
end