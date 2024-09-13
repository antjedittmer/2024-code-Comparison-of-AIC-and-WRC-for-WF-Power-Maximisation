function d2f_db2 = second_derivative_b(a1, a2, a3, c12, c13, cp3)

    % If cp3 is not provided, calculate it as 4*a3*(1 - a3)^2
    if nargin < 6
        cp3 = 4 * a3 * (1 - a3)^2;
    end

    % Compute the square root term f_sqrt(a_2)
    f_sqrt_a2 = sqrt(a1^2 * c13^2 + a2^2 * c12^2);
    
    % First term: 4(-4 + 6*a2) * (1 - 2*a1*c12)^3
    term1 = 4 * (-4 + 6 * a2) * (1 - 2 * a1 * c12)^3;
    
    % Second term: -6 * cp3 * c12^2 * (1 - 2 * f_sqrt_a2)^2 * f_sqrt_a2^(-1)
    term2 = -6 * cp3 * c12^2 * (1 - 2 * f_sqrt_a2)^2 / f_sqrt_a2;
    
    % Third term: (24 * a2^2 * cp3 * c12^4 * (1 - 2 * f_sqrt_a2)) / f_sqrt_a2^3
    term3 = (24 * a2^2 * cp3 * c12^4 * (1 - 2 * f_sqrt_a2)) / f_sqrt_a2^3;
    
    % Combine the terms to get the second derivative with respect to b
    d2f_db2 = term1 + term2 + term3;
end
