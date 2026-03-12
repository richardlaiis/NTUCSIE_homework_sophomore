function L = chol_twolevel(A)
    n = size(A, 1);
    L = zeros(n);
    for k = 1:n
        alpha = A(k, k);
        L(k, k) = sqrt(alpha);
        if k < n
            v = A(k+1:n, k);
            L(k+1:n, k) = v / L(k, k);
            for i = 1:n-k
                for j = 1:n-k
                    A(k+i, k+j) = A(k+i, k+j) - (v(i) * v(j)) / alpha;
                end
            end
        end
    end
end