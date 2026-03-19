% Define the test sizes
n_values = [10, 100, 1000, 2000, 4000];

for n = n_values
    fprintf('--- Testing for n = %d ---\n', n);
    A = generate_spd(n);
    
    % Validation
    L_builtin = chol(A, 'lower'); 
    L_vec = chol_vectorized(A);
    error_norm = norm(L_builtin - L_vec, 'fro');
    fprintf('Validation Error: %e\n', error_norm);
    
    
    tic;
    chol_vectorized(A);
    t_vec = toc;
    fprintf('Vectorized time: %f seconds\n', t_vec);
    
    tic; 
    chol_oneloop(A); 
    t_one = toc;
    fprintf('One-loop time:   %f seconds\n', t_one);
    
    tic; 
    chol_twoloop(A); 
    t_two = toc;
    fprintf('Two-loop time:   %f seconds\n', t_two);
    
    fprintf('\n');
end
