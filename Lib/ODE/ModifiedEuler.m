function w = ModifiedEuler(f,w0,a,b,n,L)
    h = (b-a)/n;
    t = @(i) a + h*(i-1);
    if nargin==5
        w = nan(1,n); w(1) = w0;
        for i = 1:n-1
            w(i+1) = w(i) + h*f(t(i)/2, w(i) + h*f(t(i),w(i))/2);
        end
    else
        w = nan(L,n); w(:,1) = w0;
        for i = 1:n-1
            w(:,i+1) = w(:,i) + h*f(t(i)/2, w(:,i) + h*f(t(i),w(:,i))/2);
        end
    end
end