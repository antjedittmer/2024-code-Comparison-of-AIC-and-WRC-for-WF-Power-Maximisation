function d2f_da3_sq = second_partial_a3(a1, a2, a3, c12, c13)
    % Calculate the square root term
    sqrt_term = sqrt(a1^2 * c13^2 + a2^2 * c12^2);
    
    % Calculate the factor (1 - 2 * sqrt_term)^3
    factor = (1 - 2 * sqrt_term)^3;
    
    % Compute the second partial derivative
    d2f_da3_sq = 4 * (-4 + 6 * a3) * factor;
end