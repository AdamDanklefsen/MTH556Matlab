function x = Solve5diag(AA,b,m,n)
    X = m+1; Y = n+1;
    U = nan(X,X,Y);
    L = nan(X,X,Y);
    x = nan(size(b));
    y = nan(size(b));
    
    U(:,:,1) = A(1);
    y(1:X) = b(1:X);
    
    for i = 2:Y
       L(:,:,i) = C(i)*U(:,:,i-1)^-1;
       U(:,:,i) = A(i) - L(:,:,i)*B(i-1);
       y((1:X) + X*(i-1)) = b((1:X) + X*(i-1)) - L(:,:,i)*y((1:X) + X*(i-2));
    end
    x((1:X) + X*(Y-1)) = U(:,:,Y)^-1*y((1:X) + X*(Y-1));
    for i = Y-1:-1:1
        x((1:X) + X*(i-1)) = U(:,:,i)^-1*(y((1:X) + X*(i-1)) - B(i)*x((1:X) + X*i));
    end
    
    
    function aa = A(i)
        aa = diag(AA((1:X) + X*(i-1),3)) + ...
            diag(AA((2:X) + X*(i-1),2),-1) + ...
            diag(AA((1:(X-1)) + X*(i-1),4),1);
    end
    function bb = B(i)
        id = (1:X) + X*(i-1);
        bb = diag(AA(id,5));
    end
    function cc = C(i)
        id = (1:X) + X*(i-1);
        cc = diag(AA(id,1));
    end
end

