function L = chol_onelevel(A)
    n = size(A, 1);
    L = zeros(n);
    for k = 1:n
        alpha = A(k, k);
        L(k, k) = sqrt(alpha);
        if k < n
            v = A(k+1:n, k);
            L(k+1:n,k) = v / L(k, k);
            for i = 1:n-k
                A(k+1:n, k+i) = A(k+1:n, k+i) - v * (v(i) / alpha);
            end
        end
    end
end