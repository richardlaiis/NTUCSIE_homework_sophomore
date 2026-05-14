% quad_least_square.m
function [a, b, c] = quad_least_square(x, y)
    x = x(:);
    y = y(:);
    
    m = length(x);

    sum_x  = sum(x);
    sum_x2 = sum(x.^2);
    sum_x3 = sum(x.^3);
    sum_x4 = sum(x.^4);
    
    sum_y   = sum(y);
    sum_xy  = sum(x .* y);
    sum_x2y = sum((x.^2) .* y);
    
    A = [sum_x4, sum_x3, sum_x2;
         sum_x3, sum_x2, sum_x;
         sum_x2, sum_x,  m];
    z = [sum_x2y; 
         sum_xy; 
         sum_y];
         
    coeffs = A \ z;
    
    a = coeffs(1);
    b = coeffs(2);
    c = coeffs(3);
end