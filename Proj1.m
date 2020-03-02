clear; close all; clc;
addpath(genpath('./Lib/ODE/'));

f  = @(t,y) -y;
fp = @(t,y) -1;
n = 2^10 + 1;
t = linspace(0,1,n);
w0 = ForwardEuler(f,.5,0,1,n);
w1 = ModifiedEuler(f,.5,0,1,n);
w2 = RK2(f,.5,0,1,n);
w3 = RK4(f,.5,0,1,n);
w4 = CrankNicholson(f,fp,.5,0,1,n);

figure(1); hold on;
plot(t,exp(-t)/2);
plot(t,w4);

f = @(t,y) -5*y + 5*t^2 + 2*t;
n = 10 + 1;
t = linspace(0,1,n);
w0 = ForwardEuler(f,1/3,0,1,n);
w1 = ModifiedEuler(f,1/3,0,1,n);
w2 = RK2(f,1/3,0,1,n);
w3 = RK4(f,1/3,0,1,n);

%figure(2); hold on;
%plot(t,exp(-5*t)/3 + t.^2);
%plot(t,w0);


f = @(t,y) -sin(y(1));
n = 2^10 + 1;
t = linspace(0,4*pi,n);
f = @(t,y) [y(2); -sin(y(1))];
w0 = RK3(f,[.5;1.2],0,4*pi,n,2);
w1 = ForwardEuler(f,[.5;1.2],0,4*pi,n,2);

figure(3); hold on;
plot(t,w0(1,:));
plot(t,w1(1,:));
plot(t,2*asin(0.649006*jacobiSN(0.395289 + t, 0.421209)));



f = @(t,y) -y;
n = 50 + 1;
t = linspace(0,1,n);
w0 = AB3(f,.5,0,1,n);
w1 = RK4(f,.5,0,1,n);
figure(4); hold on;
plot(t,w0);
%plot(t,w1);
plot(t,exp(-t)/2);

a=1;b=1;r=20;s=10;
f = @(t,y) [(r*y(1) - a*y(1)*y(2)); (-s*y(2) + a*b*y(1)*y(2))];
n = 301;
t = linspace(0,3,n);
w0 = AB3(f,[10;10],0,3,n,2);
figure(5); hold on;
plot(w0(1,:),w0(2,:));
scatter(s/a/b,r/a);
