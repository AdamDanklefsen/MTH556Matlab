function w = CrankNicholson(f,w0,a,b,h)
    n = length(a:h:b);
    w = nan(1,n); w(1) = w0;
    t = @(i) a + h*(i-1);
    for i = 1:n-1
        g  = @(x) w(i) - x + (h/2)*( f(t(i),w(i)) + f(t(i+1),x) );
        w(i+1) = fzero(g,w(i));
    end
end