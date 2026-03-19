function L = chol_oneloop(A)
    n = size(A, 1);
    L = zeros(n);
    for k = 1:n
        alpha = A(k,k);
        L(k,k) = sqrt(alpha);
        if k < n
            v = A(k+1:n, k);
            L(k+1:n, k) = v / L(k,k);
            % One-level loop: Update each column j
            for j = k+1:n
                A(k+1:n, j) = A(k+1:n, j) - v * (v(j-k) / alpha);
            end
        end
    end
end