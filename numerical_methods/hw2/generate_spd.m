function A = generate_spd(n)
    % haha 
    X = randn(n); 
    A = X * X' + n * eye(n);
endfunction

