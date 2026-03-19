function L = chol_vectorized(A)
    n = size(A, 1);
    L = zeros(n);
    for k = 1:n
        alpha = A(k,k);
        L(k,k) = sqrt(alpha);
        if k < n
            v = A(k+1:n, k);
            L(k+1:n, k) = v / L(k,k);
            % Vectorized update: No inner loops
            A(k+1:n, k+1:n) = A(k+1:n, k+1:n) - (v * v') / alpha;
        end
    end
end