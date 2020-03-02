function w = AB3(f,w0,a,b,n,L)
    h = (b-a)/n;
    t = @(i) a + h*(i-1);
    if nargin==5
        w = nan(1,n); w(1) = w0;
        w(1:3) = RK4(f,w0,a,t(4),3);
        for i = 3:n-1
            w(i+1) = w(i) + h*( 23*f(t(i),w(i)) + -16*f(t(i-1),w(i-1)) + 5*f(t(i-2),w(i-1)) )/12;
        end
    else
        w = nan(L,n); w(:,1) = w0;
        w(:,1:3) = RK4(f,w0,a,t(4),3,L);
        for i = 3:n-1
            w(:,i+1) = w(:,i) + h*( 23*f(t(i),w(:,i)) + -16*f(t(i-1),w(:,i-1)) + 5*f(t(i-2),w(:,i-1)) )/12;
        end
    end
end