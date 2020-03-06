clear; close all; clc;
addpath(genpath('./Lib/ODE/'));
%% Problem 1
f  = @(t,y) -y;
sol = @(t) exp(-t)/2;
m = 3:11;
h = 2.^-m;
t1=0;t2=2;
E = nan(length(m),3);
for i = 1:length(h)
    w1 = ModifiedEuler(f,.5,t1,t2,h(i));
    w2 = RK2(f,.5,t1,t2,h(i));
    w3 = RK4(f,.5,t1,t2,h(i));
    E(i,:) = abs(sol(2) - [w1(end);w2(end);w3(end)]);
end
a = log2(E(1:end-1,:)./E(2:end,:));
h = h(1:end-1);
A = [h' E(1:end-1,1) a(:,1) E(1:end-1,2) a(:,2) E(1:end-1,3) a(:,3)];
writematrix(A,'P1.csv');

figure(1);
subplot(2,2,1);
loglog(h,E(1:end-1,1)); hold on; loglog(h,h.^2); title('Modified Euler O(h^2)');
subplot(2,2,2);
loglog(h,E(1:end-1,2)); hold on; loglog(h,h.^2); title('RK2 O(h^2)');
subplot(2,2,3);
loglog(h,E(1:end-1,3)); hold on; loglog(h,h.^4); title('RK4 O(h^4)');


f = @(t,y) -5*y + 5*t^2 + 2*t;
sol = @(t) t.^2 + exp(-5*t)/3;
h = .1; t = 0:h:1;
w1 = ModifiedEuler(f,1/3,0,1,h);
w2 = RK2(f,1/3,0,1,h);
w3 = RK4(f,1/3,0,1,h);
figure(2);
subplot(2,2,1); hold on;
plot(t,sol(t)); plot(t,w1); title('Modified Euler');
subplot(2,2,2); hold on;
plot(t,sol(t)); plot(t,w2); title('RK2');
subplot(2,2,3); hold on;
plot(t,sol(t)); plot(t,w3); title('RK4');



%% Problem 2
f  = @(t,y) -y;
sol = @(t) exp(-t)/2;
m = 3:11;
h = 2.^-m;
t1=0;t2=2;
E = nan(length(m),1);
for i = 1:length(m)
    w0 = CrankNicholson(f,.5,t1,t2,h(i));
    E(i) = abs(sol(2) - w0(end));
end
a = log2(E(1:end-1)./E(2:end));
h = h(1:end-1);
figure(3); loglog(h,E(1:end-1,1)); hold on; loglog(h,h.^2); title('O(h^2)');
A = [h' E(1:end-1) a];
writematrix(A,'P2.csv');

f = @(t,y) -y^2 - 5*exp(-t)*sin(5*t);
h = .02; t = 0:h:2;
w0 = CrankNicholson(f,1,0,2,h);
figure(4); hold on;
plot(t,w0);

%% Problem 3
f  = @(t,y) -y;
sol = @(t) exp(-t)/2;
m = 3:11;
h = 2.^-m;
t1=0;t2=2;
E = nan(length(m),1);
for i = 1:length(m)
    w0 = RK3(f,.5,t1,t2,h(i));
    E(i) = abs(sol(2) - w0(end));
end
a = log2(E(1:end-1)./E(2:end));
h = h(1:end-1);
figure(5); loglog(h,E(1:end-1,1)); hold on; loglog(h,h.^3); title('O(h^3)');
A = [h' E(1:end-1) a];
writematrix(A,'P3.csv');

f = @(t,y) [y(2); -sin(y(1))];
h = .1; t = 0:h:10;
w0 = RK3(f,[.5;1.2],0,10,h,2);
figure(6); hold on;
plot(t,w0(1,:)); plot(t,w0(2,:));
legend('Angle','Velocity');

%% Problem 4
f  = @(t,y) -y;
sol = @(t) exp(-t)/2;
m = 3:11;
h = 2.^-m;
t1=0;t2=2;
E = nan(length(m),1);
for i = 1:length(m)
    w0 = AB3(f,.5,t1,t2,h(i));
    E(i) = abs(sol(2) - w0(end));
end
a = log2(E(1:end-1)./E(2:end));
h = h(1:end-1);
figure(7); loglog(h,E(1:end-1,1)); hold on; loglog(h,h.^3); title('O(h^3)');
A = [h' E(1:end-1) a];
writematrix(A,'P4.csv');


a=1;b=1;r=20;s=10;
f = @(t,y) [(r*y(1) - a*y(1)*y(2)); (-s*y(2) + a*b*y(1)*y(2))];
h = .01;
t = 0:h:3;
w0 = AB3(f,[10;10],0,3,h,2);
[X,Y] = meshgrid(0:1:25, 0:2:50);
Vx = r*X-a*X.*Y; Vy = -s*Y + a*b*X.*Y;
figure(12); hold on;
plot(w0(1,:),w0(2,:));
quiver(X,Y,Vx,Vy);
scatter(s/a/b,r/a);
axis([0 25 0 50]);
