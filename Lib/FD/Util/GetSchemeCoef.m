function [coef,order] = GetSchemeCoef(g,L)
    A = g.^0;
    for i = 1:length(g)-1
        A = [A;g.^i/factorial(i)];
    end
    b = zeros(size(g')); b(L+1) = 1;
    coef = A\b;
    i = length(g); order = length(g)-1;
    while(1)
        if(dot(g.^i,coef)==0)
            order = order+1;
        else
            break;
        end
    end
end

