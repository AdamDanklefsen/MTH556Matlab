%%%  HW2
%
% MTH 556
% 
% Adam Danklefsen
clear; close all; clc;
CodeCommentStart = nan; CodeCommentEnd = nan;
%% 1
CodeCommentStart;
A = [2 1; 1 4];
% a
[D,L,U] = get_jac(A);
T_jac = D^-1*(L+U);
T_gs = (D-L)^-1*U;
w = 1.5;
T_sor = (D/w-L)^-1*((1/w-1)*D+U);
disp('T_Jac = '); disp(T_jac);
disp('T_GS = ');  disp(T_gs);
disp('T_SOR = '); disp(T_sor);
% b
rho_jac = spect(T_jac);
rho_gs  = spect(T_gs );
fprintf('Jacobi Method will converge: %s\n',string(rho_jac < 1));
fprintf('Gauss-Seidel Method will converge: %s\n',string(rho_gs < 1));
fprintf('Both spectral radii are less than 1, so both methods will converge\n\n');
% c
ww = .01:.01:2;
rho = nan(size(ww));
for i = 1:length(ww)
    rho(i) = spect((D/ww(i)-L)^-1*((1/ww(i)-1)*D+U));
end
figure(1); hold on;
plot(ww,rho);
[~,i] = min(rho);
scatter(ww(i),rho(i));
text(1.2,.1,'$$(\omega = 1.04, \rho = 0.04)$$','Interpreter','latex');
xlabel('\omega'); ylabel('\rho(T_{SOR})');
title('\rho(T_{SOR}) vs \omega');

% d
w_opt = 2*(1 + sqrt(1 - spect(T_jac)^2))^-1;
fprintf('w_opt = %.4f is very close to the previously calculated value of w = %.4f\n',w_opt,ww(i));
CodeCommentEnd;
%% 2
A = [4 -1 0; -1 4 -1; 0 -1 4];
b = [2 4 10]';
e = 1e-7;
x0 = zeros(size(b));
% jacobi
tic; 
[D,L,U] = get_jac(A);
T = zeros(size(A));
c = zeros(size(b));
for i = 1:length(A)
    for j = 1:length(A)
        for k = 1:length(A)
            T(i,j) = T(i,j) + (L(k,j) + U(k,j)) * invdiag(i,k,D);
        end
    end
end
for i = 1:length(b)
    for k = 1:length(b)
        c(i) = c(i) + b(k) * invdiag(i,k,D);
    end
end
t_a = toc; tic;

err = 1; n = 0;
xnew = x0; xold = x0;
errvec = err;
while( err > e)
    n = n+1;
    for i = 1:length(T)
        for j = 1:length(T)
            xnew(i) = xnew(i) + T(i,j)*xold(j) + c(j)*(i==j);
        end
    end    
    err = norm(xnew - xold);
    xold = xnew; xnew = zeros(size(xnew));
    errvec = [errvec err];
end
t_a(2) = toc;
err_a = errvec;

% GS
tic; 
[D,L,U] = get_jac(A);
T = zeros(size(A));
B = (D-L)^-1;
c = zeros(size(b));
for i = 1:length(A)
    for j = 1:length(A)
        for k = 1:length(A)
            T(i,j) = T(i,j) + B(i,k)*U(k,j);
        end
    end
end
for i = 1:length(b)
    for k = 1:length(b)
        c(i) = c(i) + B(i,k)*b(k);
    end
end
t_b = toc; tic;

err = 1; n = 0;
xnew = x0; xold = x0;
errvec = err;
while( err > e)
    n = n+1;
    for i = 1:length(T)
        for j = 1:length(T)
            xnew(i) = xnew(i) + T(i,j)*(xold(j)*(j>=i) + xnew(j)*(i>j)) + c(j)*(i==j);
        end
    end
    err = norm(xnew - xold);
    xnew'
    xold = xnew; xnew = zeros(size(xnew));
    errvec = [errvec err];
end
t_b(2) = toc;
err_b = errvec;


% SOR
tic; 
[D,L,U] = get_jac(A);
w = 2*(1 + sqrt(7/8))^-1;
T = zeros(size(A));
B = (D/w-L)^-1;
c = zeros(size(b));
for i = 1:length(A)
    for j = 1:length(A)
        for k = 1:length(A)
            T(i,j) = T(i,j) + B(i,k) * ( (1/w-1)*D(k,j) + U(k,j) );
        end
    end
end
for i = 1:length(b)
    for k = 1:length(b)
        c(i) = c(i) + B(i,k)*b(k);
    end
end
t_c = toc; tic;

err = 1; n = 0;
xnew = x0; xold = x0;
errvec = err;
while( err > e)
    n = n+1;
    for i = 1:length(T)
        for j = 1:length(T)
            xnew(i) = xnew(i) + T(i,j)*(xold(j)*(j>=i) + xnew(j)*(i>j)) + c(i)*(i==j);
        end
    end
    err = norm(xnew - xold);
    xold = xnew; xnew = zeros(size(xnew));
    errvec = [errvec err];
end
t_c(2) = toc;
err_c = errvec;

% d
CodeCommentStart;
figure(2); hold on;
plot(1:length(err_a),err_a);
plot(1:length(err_b),err_b);
plot(1:length(err_c),err_c);
set(gca,'YScale','log');
plot([1 20], [1 1]*1e-7);
xlabel('Number of iterations');
ylabel('|| x^{new} - x^{old} ||');
title('Error Plot');
legend('Jacobi Method','Gauss-Seidel Method','SOR Method','\epsilon = 10^{-7}');
CodeCommentEnd;
function d = invdiag(a,b,D)
    if(a==b) d = 1/D(a,b); return; end
    d = 0;
end
function rho = spect(A)
    rho = max(abs(eig(A)));
end
function [D,L,U] = get_jac(A)
    D = diag(diag(A));
    L = -tril(A-D);
    U = D-A-L;
end
function a = isDiagDom(A)
    if min(size(A)) ~= length(A) a = false; return; end
    a = true;
    for i = 1:length(A)
        a = a && (abs(A(i,i)) >= sum(abs(A(i,:))) - abs(A(i,i)));
    end

end