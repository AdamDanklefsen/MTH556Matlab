function alpha = getAlpha(A,b,d,x)
    r = A * x - b;
	alpha = (r' * A * d) / (d' * A * d);
end

