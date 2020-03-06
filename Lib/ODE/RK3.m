function w = RK3(f,w0,a,b,h,L)
    n = length(a:h:b);
    t = @(i) a + h*(i-1);
    if nargin==5
    	w = nan(1,n); w(1) = w0;
        for i = 1:n-1
            k1 = f(t(i)      , w(i)              );
            k2 = f(t(i) + h/2, w(i) + h*k1/2     );
            k3 = f(t(i) + h  , w(i) - h*k1+2*h*k2);
            w(i+1) = w(i) + h*(k1 + 4*k2 + k3)/6;
        end
    else
        w = nan(L,n); w(:,1) = w0;
        for i = 1:n-1
            k1 = f(t(i)      , w(:,i)              );
            k2 = f(t(i) + h/2, w(:,i) + h*k1/2     );
            k3 = f(t(i) + h  , w(:,i) - h*k1+2*h*k2);
            w(:,i+1) = w(:,i) + h*(k1 + 4*k2 + k3)/6;
        end
    end
end