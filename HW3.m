clear; close all; clc;
addpath(genpath('./Lib/Linear/'));
addpath(genpath('./Lib/NonLinear/'));
addpath(genpath('./Lib/FD/Util/'));
A = [4 -1 0; -1 4 -1; 0 -1 4];
b = [2 4 10]';
x0 = zeros(size(b));

%% Conjugate Gradient
[x, errvec] = ConjGradMethod(A,b,x0);        % ./Lib/Linear/
errvec = errvec(errvec~=0); errvec = errvec(errvec~=1);
fprintf('Conjugate Gradient Method Converged to x = [ '); fprintf(repmat('%f ',1,length(x)),x);
fprintf('] in %d steps\n',length(errvec));

%% Newton's Method
f = @(x) exp(-x)-x; fp = @(x) -exp(-x)-1;
[x, errvec] = NewtonsMethod(f,fp,5.25,100);   % ./Lib/NonLinear/
errvec = errvec(errvec~=0); errvec = errvec(errvec~=1);
fprintf('Newton''s Method Converged to x = %.20f in %d steps\n',x,length(errvec));


%% Error Analysis
f = @(x) exp(x);
dt = .001;
dt = dt:dt:.1; x = 1;
dydt = FD_012_1(f,x,dt);
figure(1);
loglog(dt,abs(dydt-exp(1))); hold on;
loglog(dt,5*dt.^2);
xlabel("$h$",'Interpreter','latex');
ylabel("$|Err|$",'Interpreter','latex');
title('Error Plot');
legend('Error','Cdt^2'); legend('Location','Northwest');


function dydt = FD_012_1(f,x,dt)
    a = 0:2; L = 1;
    g = GetSchemeCoef(a,L);     % ./Lib/FD/Util/
    dydt = 0;
    for i = 1:length(a)
        dydt = dydt + g(i)*f(x+a(i)*dt);
    end
    dydt = dydt.*dt.^-L;
end