function [D,L,U] = getDLU(A)
    D = diag(diag(A));
    L = -tril(A-D);
    U = D-A-L;
end

