function L = chol_twoloop(A)
    n = size(A, 1);
    L = zeros(n);
    for k = 1:n
        alpha = A(k,k);
        L(k,k) = sqrt(alpha);
        if k < n
            v = A(k+1:n, k);
            L(k+1:n, k) = v / L(k,k);
            % Two-level loop: Update each element B(i,j)
            for j = k+1:n
                for i = k+1:n
                    A(i, j) = A(i, j) - v(i-k) * (v(j-k) / alpha);
                end
            end
        end
    end
end