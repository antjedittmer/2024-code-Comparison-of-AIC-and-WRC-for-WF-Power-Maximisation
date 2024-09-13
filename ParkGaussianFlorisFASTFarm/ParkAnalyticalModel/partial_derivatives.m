function [df_da1, df_da2, df_da3] = partial_derivatives(a1, a2, a3, c12, c13)
    % Partial derivative with respect to a1
    term1_a1 = 4 * (1 - 4 * a1 + 3 * a1^2);
    term2_a1 = -24 * a2 * (1 - a2)^2 * c12 * (1 - 2 * a1 * c12)^2;
    sqrt_term_a1 = sqrt(a1^2 * c13^2 + a2^2 * c12^2);
    term3_a1 = -24 * a1 * a3 * (1 - a3)^2 * (c13^2 * (1 - 2 * sqrt_term_a1)^2) / sqrt_term_a1;
    df_da1 = term1_a1 + term2_a1 + term3_a1;

    % Partial derivative with respect to a2
    term1_a2 = 4 * (1 - 4 * a2 + 3 * a2^2) * (1 - 2 * a1 * c12)^3;
    sqrt_term_a2 = sqrt(a1^2 * c13^2 + a2^2 * c12^2);
    term2_a2 = -24 * a2 * a3 * (1 - a3)^2 * (c12^2 * (1 - 2 * sqrt_term_a2)^2) / sqrt_term_a2;
    df_da2 = term1_a2 + term2_a2;

    % Partial derivative with respect to a3
    sqrt_term_a3 = sqrt(a1^2 * c13^2 + a2^2 * c12^2);
    df_da3 = 4 * (1 - 4 * a3 + 3 * a3^2) * (1 - 2 * sqrt_term_a3)^3;
end
