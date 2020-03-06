function w = AB3(f,w0,a,b,h,L)
    n = length(a:h:b);
    t = @(i) a + h*(i-1);
    if nargin==5
        w = nan(1,n); w(1) = w0;
        w(1:3) = RK4(f,w0,a,t(3),h);
        for i = 3:n-1
            w(i+1) = w(i) + h*( 23*f(t(i),w(i)) + -16*f(t(i-1),w(i-1)) + 5*f(t(i-2),w(i-2)) )/12;
        end
    else
        w = nan(L,n); w(:,1) = w0;
        w(:,1:3) = RK4(f,w0,a,t(3),h,L);
        for i = 3:n-1
            w(:,i+1) = w(:,i) + h*( 23*f(t(i),w(:,i)) + -16*f(t(i-1),w(:,i-1)) + 5*f(t(i-2),w(:,i-2)) )/12;
        end
    end
end