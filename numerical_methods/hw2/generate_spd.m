function A = generate_spd(n)
    X = randn(n); 
    A = X * X' + n * eye(n);
endfunction

