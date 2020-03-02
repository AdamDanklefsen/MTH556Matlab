% [Vector, vector] = JacobiMethod(Matrix, Vector, Vector)
% Runs FixedPointLinear with Jacobi Paramters
% @ Matrix A
% @ Vector b
% @ Vector x0
% 
% @ Vector X            - Solution
% @ vector errvec       - array of error values

function [X,errvec] = JacobiMethod(A,b,x0)
    addpath(genpath('./Utility/'));
    [D,L,U] = getDLU(A);
    T = D^-1*(L+U);
    c = D^-1*b;
    [X,errvec] = FixedPointLinear(T,c,x0);
end

