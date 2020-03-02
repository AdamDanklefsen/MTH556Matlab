function [x, errvec] = NewtonsMethod(f,fp,x0,Nmax)
    tol = 1e-14;
    err = 1; errvec = err;
    k = 0; x = x0;
    while(err>tol && k<Nmax)
        k=k+1;
        x = x - f(x)/fp(x);
        err = f(x);
        errvec = [errvec err];
    end
end