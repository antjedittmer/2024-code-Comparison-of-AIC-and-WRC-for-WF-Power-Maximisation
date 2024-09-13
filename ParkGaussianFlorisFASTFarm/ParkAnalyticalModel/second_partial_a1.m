function d2f_da1_sq = second_partial_a1(a1, a2, a3, c12, c13)
    % First term: 4(-4 + 6a1)
    first_term = 4 * (-4 + 6 * a1);
    
    % Second term: differentiate with respect to a1 (chain rule)
    second_term = -24 * a2 * (1 - a2)^2 * c12 * (2 * (1 - 2 * a1 * c12) * (-2 * c12));
    
    % Third term: complex part with the square root and product rule
    sqrt_term = sqrt(a1^2 * c13^2 + a2^2 * c12^2);
    factor = (1 - 2 * sqrt_term)^2;
    third_term = -24 * a3 * (1 - a3)^2 * (c13^2 * factor) / sqrt_term;
    
    % Product rule for third term, derivative with respect to a1
    derivative_sqrt = 24 * a1 * a3 * (1 - a3)^2 * c13^2 * (1 - 2 * sqrt_term)^2 / (sqrt_term^2);
    
    % Combine all terms
    d2f_da1_sq = first_term + second_term + third_term + derivative_sqrt;
end
