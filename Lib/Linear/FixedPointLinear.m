% [Vector, vector] = FixedPointLinear(Matrix,Vector,Vector,bool,double)
% Solves Fixed Point problem x_k+1 = Tx_k+c
% @ Matrix T
% @ Vector b
% @ Vector x0           - initial x value
% @ opt bool useNewX    - use previously calculated x values
% @ opt double tol      - tolerance
% 
% @ Vector ans          - Solution
% @ vector errvec       - array of error values

function [X, errvec] = FixedPointLinear(T,c,x0,useNewX,tol)
    switch(nargin)
        case 4
            if( ~(useNewX==0 || useNewX==1) )
                tol = useNewX;
                useNewX = false;
            else
                tol = 1e-7;
            end
        case 3
            useNewX = false;
            tol = 1e-7;
    end
    
    Xnew = x0; Xold = x0; err = 1; errvec = err;
    if(useNewX)
        while(err > tol)
            for i = 1:length(T)
                a = T*([Xnew(1:i-1)' Xold(i:end)']') + c;
                Xnew(i) = a(i);
            end
            err = norm(Xnew-Xold);
            Xold = Xnew;
            errvec = [errvec err];
        end 
    else
        while(err > tol)
            Xnew = T*Xold + c;
            err = norm(Xnew-Xold);
            Xold = Xnew;
            errvec = [errvec err];
        end
    end
    
    X = Xnew;
end