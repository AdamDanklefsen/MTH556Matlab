function lambda = GetLambda(A,b,d,x)
    r = A * x - b;
    lambda = -(d' * r) / (d' * A * d);
end

