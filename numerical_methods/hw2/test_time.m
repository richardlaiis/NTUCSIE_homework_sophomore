n_values = [10, 100, 1000, 4000]; % Suggested n values 

for n = n_values
    A = generate_spd(n); % Create the matrix [cite: 49]
    
    % Measure Vectorized Version
    tic;
    chol_vectorized(A);
    time_vec = toc;
    
    % Measure One-level Loop Version
    tic;
    chol_onelevel(A);
    time_one = toc;

    % Measure Two-level Loop Version
    tic;
    chol_twolevel(A);
    time_two = toc;
    
    fprintf('For n=%d: Vectorized = %f s, One-loop = %f s, Two-loop = %f s\n', n, time_vec, time_one, time_two);
end