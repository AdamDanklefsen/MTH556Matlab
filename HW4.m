clear; close all; clc;
%% P1
%% P2
f = @(x) x.^2 + exp(x);
S1 = @(f,h) (2*f(1+h) + 3*f(1) - 6*f(1-h) + f(1-2*h))/6./h;
S2 = @(f,h) (f(1+2*h) - 2*f(1+h) + 2*f(1-h) - f(1-2*h))/2./h./h./h;
S3 = @(f,h) (-f(1+2*h) + 16*f(1+h) - 30*f(1) + 16*f(1-h) - f(1-2*h))/12./h./h;

e = exp(1);
f0 = 1+e;
f1 = 2+e;
f2 = 2+e;
f3 = e;
h = 2.^-(1:10);

figure(1); subplot(2,2,1);
loglog(h,abs(S1(f,h)-f1));
hold on;
loglog(h,h.^3);
title('S1 Error Order 3');

subplot(2,2,2);
loglog(h,abs(S2(f,h)-f3));
hold on;
loglog(h,h.^2);
title('S2 Error Order 2');

subplot(2,2,3);
loglog(h,abs(S3(f,h)-f2));
hold on;
loglog(h,h.^4);
title('S3 Error Order 4');

%% P3
%% P4
clear; close all; clc;
addpath(genpath('./Lib/ODE/'));

% 1
f1 = @(t,y) y+2*t*exp(2*t);
sol1 = @(t) 3*exp(t) -2*exp(2*t) + 2*t.*exp(2*t);
f2 = @(t,y) -y+2*t*exp(-t);
sol2 = @(t) t.^2.*exp(-t) + exp(-t);

h = .05*2.^-(1:16);
y2 = nan(2,length(h));
E  = nan(2,length(h));
for i = 1:length(h)
    t = 0:h(i):2; n = length(t);
    w1 = ForwardEuler(f1,1,0,2,n);
    w2 = ForwardEuler(f2,1,0,2,n);
    y2(:,i) = [w1(end); w2(end)];
    E(:,i)  = abs([sol1(t(end)); sol2(t(end))] - [w1(end); w2(end)]);
end

fprintf('h\t\t\ty1(2)\t\ty2(2)\n');
for i = 1:length(h)
    fprintf('.05*2^-%d\t%.6f\t%.6f\n',i,y2(1,i),y2(2,i));
    
end

figure(4); loglog(h,E(1,:)); hold on; loglog(h,h.^1); title('y_1(t) Error Plot Order 1');
figure(5); loglog(h,E(2,:)); hold on; loglog(h,h.^1); title('y_2(t) Error Plot Order 1');



