% [Vector, vector] = SORMethod(Matrix, Vector, Vector)
% Runs FixedPointLinear with Jacobi Paramters
% @ Matrix A
% @ Vector b
% @ Vector x0
% @ Scalar w
% 
% @ Vector X            - Solution
% @ vector errvec       - array of error values

function [X,errvec] = SORMethod(A,b,x0,w)
    addpath(genpath('./Utility/'));
    [D,L,U] = getDLU(A);
    T = (D*(1/w)-L)^-1*((1/w-1)*D+U);
    c = (D*(1/w)-L)^-1*b;
    [X,errvec] = FixedPointLinear(T,c,x0,true);
end

