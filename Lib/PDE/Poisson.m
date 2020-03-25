function [AA,b,nds] = Poisson(R,m,n,f,BC)
dim = (m+1)*(n+1);
x = @(i) R(1) +  (i-1)*(R(2)-R(1))/m;
y = @(j) R(3) +  (j-1)*(R(4)-R(3))/n;
dx = x(1)-x(0); dy = y(1)-y(0);
dxx = dx*dx; A = (dx/dy)^2;
BCL = BC(1); BCR = BC(2);
BCD = BC(3); BCU = BC(4);
AA = nan(dim,5);
b = nan(dim,1);
nds(dim) = "";

for j = 1:n+1
    for i = 1:m+1
        id = (j-1)*(n+1)+i;
        [row, B, nd] = genRow(i,j);
        AA(id,:) = row;
        b(id) = B;
        nds(id) = nd;
    end
end


    function [M, B, nd] = genRow(i,j)
        M = nan(1,5);
        % Center
        if(2<=i && i<=m && 2<=j && j<=n)
            M(1) = A;
            M(2) = 1;
            M(3) = -2*(1+A);
            M(4) = 1;
            M(5) = A;
            B = dxx*f(x(i),y(j));
            nd = 'Center';
            return;
        end
        % Corners
        if(i==1 && j==1) % DL
            nd = 'DL';
            if(BCL.Dirichlet)
                M = [0 0 1 0 0];
                B = BCL.L(x(i),y(j))/BCL.a(x(i),y(j));
                return;
            end
            if(BCD.Dirichlet)
                M = [0 0 1 0 0];
                B = BCD.L(x(i),y(j))/BCD.a(x(i),y(j));
                return;
            end
            M    = [0 0 nan 2 2*A];
            M(3) = -2*(1+A + BCL.a(x(i),y(j))*dx/BCL.b(x(i),y(j)) + A*BCD.a(x(i),y(j))*dy/BCD.b(x(i),y(j)));
            B    = dxx*f(x(i),y(j)) - BCL.L(x(i),y(j))*2*dx/BCL.b(x(i),y(j)) - A*BCD.L(x(i),y(j))*2*dy/BCD.b(x(i),y(j));
            return;
        end
        if(i==m+1 && j==1) % DR
            nd = 'DR';
            if(BCR.Dirichlet)
                M = [0 0 1 0 0];
                B = BCR.L(x(i),y(j))/BCR.a(x(i),y(j));
                return;
            end
            if(BCD.Dirichlet)
                M = [0 0 1 0 0];
                B = BCD.L(x(i),y(j))/BCD.a(x(i),y(j));
                return;
            end
            M    = [0 2 nan 0 2*A];
            M(3) = -2*(1+A + BCR.a(x(i),y(j))*dx/BCR.b(x(i),y(j)) + A*BCD.a(x(i),y(j))*dy/BCD.b(x(i),y(j)));
            B    = dxx*f(x(i),y(j)) - BCR.L(x(i),y(j))*2*dx/BCR.b(x(i),y(j)) - A*BCD.L(x(i),y(j))*2*dy/BCD.b(x(i),y(j));
            return;
        end
        if(i==1 && j==n+1) % UL
            nd = 'UL';
            if(BCL.Dirichlet)
                M = [0 0 1 0 0];
                B = BCL.L(x(i),y(j))/BCL.a(x(i),y(j));
                return;
            end
            if(BCU.Dirichlet)
                M = [0 0 1 0 0];
                B = BCU.L(x(i),y(j))/BCU.a(x(i),y(j));
                return;
            end
            M    = [2*A 0 nan 2 0];
            M(3) = -2*(1+A + BCL.a(x(i),y(j))*dx/BCL.b(x(i),y(j)) + A*BCU.a(x(i),y(j))*dy/BCU.b(x(i),y(j)));
            B    = dxx*f(x(i),y(j)) - BCL.L(x(i),y(j))*2*dx/BCL.b(x(i),y(j)) - A*BCU.L(x(i),y(j))*2*dy/BCU.b(x(i),y(j));
            return;
        end
        if(i==m+1 && j==n+1) % UR
            nd = 'UR';
            if(BCR.Dirichlet)
                M = [0 0 1 0 0];
                B = BCR.L(x(i),y(j))/BCR.a(x(i),y(j));
                return;
            end
            if(BCU.Dirichlet)
                M = [0 0 1 0 0];
                B = BCU.L(x(i),y(j))/BCU.a(x(i),y(j));
                return;
            end
            M    = [2*A 2 nan 0 0];
            M(3) = -2*(1+A + BCR.a(x(i),y(j))*dx/BCR.b(x(i),y(j)) + A*BCU.a(x(i),y(j))*dy/BCU.b(x(i),y(j)));
            B    = dxx*f(x(i),y(j)) - BCR.L(x(i),y(j))*2*dx/BCR.b(x(i),y(j)) - A*BCU.L(x(i),y(j))*2*dy/BCU.b(x(i),y(j));
            return;
        end
        % Left
        if(i==1)
            nd = 'Left';
            if(BCL.Dirichlet)
                M = [0 0 1 0 0];
                B = BCL.L(x(i),y(j))/BCL.a(x(i),y(j));
                return;
            end
            M(1) = A;
            M(2) = 0;
            M(3) = -2*(BCL.a(x(i),y(j))*dx/BCL.b(x(i),y(j)) + 1 + A);
            M(4) = 2;
            M(5) = A;
            B    = dxx*f(x(i),y(j)) - BCL.L(x(i),y(j))*2*dx/BCL.b(x(i),y(j));
            return;
        end
        % Right
        if(i==m+1)
            nd = 'Right';
            if(BCR.Dirichlet)
                M = [0 0 1 0 0];
                B = BCR.L(x(i),y(j))/BCR.a(x(i),y(j));
                return;
            end
            M(1) = A;
            M(2) = 2;
            M(3) = -2*(BCR.a(x(i),y(j))*dx/BCR.b(x(i),y(j)) + 1 + A);
            M(4) = 0;
            M(5) = A;
            B    = dxx*f(x(i),y(j)) - BCR.L(x(i),y(j))*2*dx/BCR.b(x(i),y(j));
            return;
        end
        % Down
        if(j==1)
            nd = 'Down';
            if(BCD.Dirichlet)
                M = [0 0 1 0 0];
                B = BCD.L(x(i),y(j))/BCD.a(x(i),y(j));
                return;
            end
            M(1) = 0;
            M(2) = 1;
            M(3) = -2*(A*BCD.a(x(i),y(j))*dy/BCD.b(x(i),y(j)) + 1 + A);
            M(4) = 1;
            M(5) = 2*A;
            B    = dxx*f(x(i),y(j)) - A*BCD.L(x(i),y(j))*2*dy/BCD.b(x(i),y(j));
            return;
        end
        % Up
        if(j==n+1)
            nd = 'Up';
            if(BCU.Dirichlet)
                M = [0 0 1 0 0];
                B = BCU.L(x(i),y(j))/BCU.a(x(i),y(j));
                return;
            end
            M(1) = 2*A;
            M(2) = 1;
            M(3) = -2*(A*BCU.a(x(i),y(j))*dy/BCU.b(x(i),y(j)) + 1 + A);
            M(4) = 1;
            M(5) = 0;
            B    = dxx*f(x(i),y(j)) - A*BCU.L(x(i),y(j))*2*dy/BCU.b(x(i),y(j));
            return;
        end
    end
end