% [Vector, vector] = Gauss_SeidelMethod(Matrix, Vector, Vector)
% Runs FixedPointLinear with Gauss-Seidel Paramters
% @ Matrix A
% @ Vector b
% @ Vector x0
% 
% @ Vector X            - Solution
% @ vector errvec       - array of error values

function [X,errvec] = Gauss_SeidelMethod(A,b,x0)
    addpath(genpath('./Utility/'));
    [D,L,U] = getDLU(A);
    T = (D-L)^-1*U;
    c = (D-L)^-1*b;
    [X,errvec] = FixedPointLinear(T,c,x0,true);
end

