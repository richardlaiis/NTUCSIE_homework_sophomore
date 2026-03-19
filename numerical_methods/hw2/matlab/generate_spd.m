function A = generate_spd(n)
    % Generate an n x n symmetric positive definite matrix
    X = randn(n); 
    A = X * X' + n * eye(n);
end