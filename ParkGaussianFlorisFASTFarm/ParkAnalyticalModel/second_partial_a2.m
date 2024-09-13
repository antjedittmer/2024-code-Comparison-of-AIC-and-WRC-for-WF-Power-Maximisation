function d2f_da2_sq = second_partial_a2(a1, a2, a3, c12, c13)
    % First term: 4(-4 + 6a2)(1 - 2a1 * c12)^3
    first_term = 4 * (-4 + 6 * a2) * (1 - 2 * a1 * c12)^3;
    
    % Second term: complex part with the square root and product rule
    sqrt_term = sqrt(a1^2 * c13^2 + a2^2 * c12^2);
    factor = (1 - 2 * sqrt_term)^2;
    second_term = -24 * a3 * (1 - a3)^2 * (c12^2 * factor) / sqrt_term;
    
    % Derivative with respect to a2 (part of product rule)
    derivative_sqrt = 24 * a2 * a3 * (1 - a3)^2 * c12^2 * (1 - 2 * sqrt_term)^2 / (sqrt_term^2);
    
    % Combine both terms
    d2f_da2_sq = first_term + second_term + derivative_sqrt;
end
