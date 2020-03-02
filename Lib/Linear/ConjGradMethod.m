function [X, errvec] = ConjGradMethod(A,b,x0)
    addpath(genpath('./Utility/'));
    Xnew = x0; Xold = x0; err = 1; errvec = err; tol = 1e-7;

	r = A * Xold - b;
	d = -r;
	Lambda = getLambda(A, b, d, Xold);
	Xnew = Xold + Lambda * d;
	alpha = nan;
	Xold = Xnew;
	err = norm(A * Xnew - b);
	errvec = [errvec err];

	while(err > tol) 
		r = A * Xold - b;
		alpha = getAlpha(A, b, d, Xold);
		d = -r + alpha * d;
		Lambda = getLambda(A, b, d, Xold);
		Xnew = Xold + Lambda * d;

		err = norm(A * Xnew - b);
		Xold = Xnew;
		errvec = [errvec err];
    end
    X = Xnew;
end

