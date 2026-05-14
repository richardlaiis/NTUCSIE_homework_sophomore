% my_spline.m
function yq = my_spline(x, y, xq)
    x = x(:);
    n = length(x) - 1; 
    
    fp_0 = y(1);         
    fp_n = y(end);       
    a = y(2:end-1);
    a = a(:);
    
    h = diff(x);
    A = zeros(n+1, n+1);
    z = zeros(n+1, 1);
    
    A(1, 1) = 2 * h(1);
    A(1, 2) = h(1);
    z(1) = (3 / h(1)) * (a(2) - a(1)) - 3 * fp_0;

    for i = 2:n
        A(i, i-1) = h(i-1);
        A(i, i)   = 2 * (h(i-1) + h(i));
        A(i, i+1) = h(i);
        
        term_right = (3 / h(i)) * (a(i+1) - a(i));
        term_left  = (3 / h(i-1)) * (a(i) - a(i-1));
        z(i) = term_right - term_left;
    end

    A(n+1, n)   = h(n);
    A(n+1, n+1) = 2 * h(n);
    z(n+1)      = 3 * fp_n - (3 / h(n)) * (a(n+1) - a(n));
    
    c = A \ z;
    
    b = zeros(n, 1);
    d = zeros(n, 1);
    for j = 1:n
        b(j) = (1 / h(j)) * (a(j+1) - a(j)) - (h(j) / 3) * (2 * c(j) + c(j+1));
        d(j) = (c(j+1) - c(j)) / (3 * h(j));
    end
    
    yq = zeros(size(xq));
    for k = 1:numel(xq)
        xq_val = xq(k);
        idx = sum(x <= xq_val);
        if idx == 0; idx = 1; elseif idx > n; idx = n; end
        dx = xq_val - x(idx);
        yq(k) = a(idx) + b(idx)*dx + c(idx)*dx^2 + d(idx)*dx^3;
    end
end