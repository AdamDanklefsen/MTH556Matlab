function w = CrankNicholson(f,fp,w0,a,b,n)
    addpath(genpath('./Lib/NonLinear/'));
    w = nan(1,n); w(1) = w0;
    h = (b-a)/n;
    t = @(i) a + h*(i-1);
    for i = 1:n-1
        g  = @(x) w(i) - x + (h/2)*( f(t(i),w(i)) + f(t(i+1),x) );
        gp = @(x) (h/2)*fp(t(i+1),x) - 1;
        w(i+1) = NewtonsMethod(g,gp,w(i),100);
    end
end