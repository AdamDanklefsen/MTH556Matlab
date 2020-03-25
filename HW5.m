clear; close all; clc;
addpath(genpath('./LIB/PDE/'));
%% 1.7
clear; close all; clc;
% setup
BCL = BC_rect(@(x,y) 1,     @(x,y) y/(1+y^2)                    );
%BCL = BC_rect(@(x,y) 0,     @(x,y) -2*y/(1+y^2)^2,              @(x,y) -1);
BCR = BC_rect(@(x,y) 1,     @(x,y) y/(4+y^2)                    );
%BCR = BC_rect(@(x,y) 0,     @(x,y) -4*y/(4+y^2)^2,              @(x,y) 1);
BCD = BC_rect(@(x,y) 1,     @(x,y) 0                            );
%BCD = BC_rect(@(x,y) 0,     @(x,y) (1+x)^-2,                    @(x,y) -1);
BCU = BC_rect(@(x,y) 1,     @(x,y) ((1+x)^2+1)^-1               );
%BCU = BC_rect(@(x,y) 0,     @(x,y) ((x+1)^2-1)*((1+x)^2+1)^-2,  @(x,y) 1);

f = @(x,y) 0;
BC = [BCL BCR BCD BCU];
R = [0 1 0 1];
SOL = @(X,Y) Y ./ ( (1+X).^2 + Y.^2 );

% loop
mm = 3:8;
E = nan(1,length(mm));
rt = nan(1,length(mm));
figure(1); hold on;
for i = 1:length(mm)
    tic;
    % step dependent setup
    MN = 2^mm(i);
    m = MN;
    n = MN;
    x = linspace(R(1),R(2),m+1);
    y = linspace(R(3),R(4),n+1);
    [X, Y] = meshgrid(x,y);
    
    % Solve
    [AA, b] = Poisson(R,m,n,f,BC);
    w = Solve5diag(AA,b,m,n);
    w = reshape(w,[m+1,n+1])';
    E(i) = norm(w-SOL(X,Y))/norm(SOL(X,Y));
    rt(i) = toc;
    disp(mm(i));
    
    if(i~=length(mm))
        subplot(2,3,i); hold on;
        imagesc(x,y,w); title(strcat('h = 2^{-', string(mm(i)), '}'));
    end
end
subplot(2,3,6); hold on;
imagesc(x,y,SOL(X,Y)); title('Exact Solution');
sgtitle('P 9.1.7');
% post
h = 2.^-mm;
a = log2(E(1:end-1)./E(2:end));
h = h(1:end-1);
E = E(1:end-1);
res = [h' E' a'];

%% 1.8
clearvars -except res; close all; clc;
BCL = BC_rect(@(x,y) 1, @(x,y) sin(pi*y));
BCR = BC_rect(@(x,y) 1, @(x,y) exp(pi)*sin(pi*y) + y^2/2);
BCD = BC_rect(@(x,y) 1, @(x,y) 0);
BCU = BC_rect(@(x,y) 1, @(x,y) x^2/2);
f = @(x,y) x^2 + y^2;
BC = [BCL BCR BCD BCU];
R = [0 1 0 1];
SOL = @(X,Y) exp(pi*X).*sin(pi*Y) + (X.*Y).^2/2;

% loop
mm = 3:8;
E = nan(1,length(mm));
rt = nan(1,length(mm));
figure(1); hold on;
for i = 1:length(mm)
    tic;
    % step dependent setup
    MN = 2^mm(i);
    m = MN;
    n = MN;
    x = linspace(R(1),R(2),m+1);
    y = linspace(R(3),R(4),n+1);
    [X, Y] = meshgrid(x,y);
    
    % Solve
    [AA, b] = Poisson(R,m,n,f,BC);
    w = Solve5diag(AA,b,m,n);
    w = reshape(w,[m+1,n+1])';
    E(i) = norm(w-SOL(X,Y))/norm(SOL(X,Y));
    rt(i) = toc;
    disp(mm(i));
    
    if(i~=length(mm))
        subplot(2,3,i); hold on;
        imagesc(x,y,w); title(strcat('h = 2^{-', string(mm(i)), '}'));
    end
end
subplot(2,3,6); hold on;
imagesc(x,y,SOL(X,Y)); title('Exact Solution');
sgtitle('P 9.1.7');
% post
h = 2.^-mm;
a = log2(E(1:end-1)./E(2:end));
h = h(1:end-1);
E = E(1:end-1);
res = [res E' a'];

%% 2.7
clearvars -except res; close all; clc;
%             alpha                 L                           beta
BCL = BC_rect(@(x,y) 1,             @(x,y) y/(1+y^2)            );
BCR = BC_rect(@(x,y) 0,             @(x,y) -4*y*(4+y^2)^-2,     @(x,y) 1);
BCD = BC_rect(@(x,y) 1,             @(x,y) 0                    );
BCU = BC_rect(@(x,y) 2/((1+x)^2+1), @(x,y) ((1+x)^2+1)^-1,      @(x,y) 1);
f = @(x,y) 0;
BC = [BCL BCR BCD BCU];
R = [0 1 0 1];
SOL = @(X,Y) Y ./ ( (1+X).^2 + Y.^2 );

% loop
mm = 3:8;
E = nan(1,length(mm));
rt = nan(1,length(mm));
figure(1); hold on;
for i = 1:length(mm)
    tic;
    % step dependent setup
    MN = 2^mm(i);
    m = MN;
    n = MN;
    x = linspace(R(1),R(2),m+1);
    y = linspace(R(3),R(4),n+1);
    [X, Y] = meshgrid(x,y);
    
    % Solve
    [AA, b] = Poisson(R,m,n,f,BC);
    w = Solve5diag(AA,b,m,n);
    w = reshape(w,[m+1,n+1])';
    E(i) = norm(w-SOL(X,Y))/norm(SOL(X,Y));
    rt(i) = toc;
    disp(mm(i));
    
    if(i~=length(mm))
        subplot(2,3,i); hold on;
        imagesc(x,y,w); title(strcat('h = 2^{-', string(mm(i)), '}'));
    end
end
subplot(2,3,6); hold on;
imagesc(x,y,SOL(X,Y)); title('Exact Solution');
sgtitle('P 9.1.7');
% post
h = 2.^-mm;
a = log2(E(1:end-1)./E(2:end));
h = h(1:end-1);
E = E(1:end-1);
res = [res E' a'];

%A = diag(AA(:,3)) + diag(AA(2:end,2),-1) + diag(AA(1:end-1,4),1) + diag(AA(10:end,1),-9) + diag(AA(1:end-9,5),9);

%% 2.9
clearvars -except res; close all; clc;
%             alpha                 L                               beta
BCL = BC_rect(@(x,y) 0,             @(x,y) -4*sin(6*y),             @(x,y) -1);
BCR = BC_rect(@(x,y) 0,             @(x,y) -4*sin(6*y+4),           @(x,y) 1);
BCD = BC_rect(@(x,y) 1,             @(x,y) cos(4*x)                 );
BCU = BC_rect(@(x,y) 1,             @(x,y) cos(4*x+6)-sin(4*x+6),   @(x,y) 1/6);
f = @(x,y) -52*cos(4*x+6*y);
BC = [BCL BCR BCD BCU];
R = [0 1 0 1];
SOL = @(X,Y) cos(4*X + 6*Y);

% loop
mm = 3:8;
E = nan(1,length(mm));
rt = nan(1,length(mm));
figure(1); hold on;
for i = 1:length(mm)
    tic;
    % step dependent setup
    MN = 2^mm(i);
    m = MN;
    n = MN;
    x = linspace(R(1),R(2),m+1);
    y = linspace(R(3),R(4),n+1);
    [X, Y] = meshgrid(x,y);
    
    % Solve
    [AA, b] = Poisson(R,m,n,f,BC);
    w = Solve5diag(AA,b,m,n);
    w = reshape(w,[m+1,n+1])';
    E(i) = norm(w-SOL(X,Y))/norm(SOL(X,Y));
    rt(i) = toc;
    disp(mm(i));
    
    if(i~=length(mm))
        subplot(2,3,i); hold on;
        imagesc(x,y,w); title(strcat('h = 2^{-', string(mm(i)), '}'));
    end
end
subplot(2,3,6); hold on;
imagesc(x,y,SOL(X,Y)); title('Exact Solution');
sgtitle('P 9.1.7');
% post
h = 2.^-mm;
a = log2(E(1:end-1)./E(2:end));
h = h(1:end-1);
E = E(1:end-1);
res = [res E' a'];

%% 2.10
clearvars -except res; close all; clc;
%             alpha                 L                               beta
BCL = BC_rect(@(x,y) 0,             @(x,y) 2/(1+y^2),               @(x,y) -1);
BCR = BC_rect(@(x,y) 1,             @(x,y) log(y^2+4)               );
BCD = BC_rect(@(x,y) 0,             @(x,y) 0,                       @(x,y) -1);
BCU = BC_rect(@(x,y) 1,             @(x,y) log(x^2+1)               );
f = @(x,y) 0;
BC = [BCL BCR BCD BCU];
R = [1 2 0 1];
SOL = @(X,Y) log(X.^2 + Y.^2);

% loop
mm = 3:8;
E = nan(1,length(mm));
rt = nan(1,length(mm));
figure(1); hold on;
for i = 1:length(mm)
    tic;
    % step dependent setup
    MN = 2^mm(i);
    m = MN;
    n = MN;
    x = linspace(R(1),R(2),m+1);
    y = linspace(R(3),R(4),n+1);
    [X, Y] = meshgrid(x,y);
    
    % Solve
    [AA, b, nds] = Poisson(R,m,n,f,BC);
    w = Solve5diag(AA,b,m,n);
    w = reshape(w,[m+1,n+1])';
    nds = reshape(nds,[m+1,n+1])';
    E(i) = norm(w-SOL(X,Y))/norm(SOL(X,Y));
    rt(i) = toc;
    disp(mm(i));
    
    if(i~=length(mm))
        subplot(2,3,i); hold on;
        imagesc(x,y,w); title(strcat('h = 2^{-', string(mm(i)), '}'));
    end
end
subplot(2,3,6); hold on;
imagesc(x,y,SOL(X,Y)); title('Exact Solution');
sgtitle('P 9.1.7');
% post
h = 2.^-mm;
a = log2(E(1:end-1)./E(2:end));
h = h(1:end-1);
E = E(1:end-1);
res = [res E' a'];



function plotBC(R,BC)
    if(BC(1).Dirichlet)
        plot([R(1) R(1)], [R(3) R(4)],'R','Linewidth',2);
    else
        plot([R(1) R(1)], [R(3) R(4)],'B','Linewidth',2);
    end
    if(BC(2).Dirichlet)
        plot([R(2) R(2)], [R(3) R(4)],'R','Linewidth',2);
    else
        plot([R(2) R(2)], [R(3) R(4)],'B','Linewidth',2);
    end
    if(BC(3).Dirichlet)
        plot([R(1) R(2)], [R(3) R(3)],'R','Linewidth',2);
    else
        plot([R(1) R(2)], [R(3) R(3)],'B','Linewidth',2);
    end
    if(BC(4).Dirichlet)
        plot([R(1) R(2)], [R(4) R(4)],'R','Linewidth',2);
    else
        plot([R(1) R(2)], [R(4) R(4)],'B','Linewidth',2);
    end
end