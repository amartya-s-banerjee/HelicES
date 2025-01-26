function [Bessel_Roots, Bessel_Evals] = Setup_Bessel_Roots_and_Evals(S)

Bessel_Roots = zeros(S.K_max, 2*S.N_max+1);
Bessel_Evals = zeros(S.K_max, 2*S.N_max+1);

% Setting up the roots b_nN^k
for n = -S.N_max:1:S.N_max
    Bessel_Roots(:, n+S.N_max+1) = besselzero_new((n*S.N), S.K_max, 1);
end

% Setting up J(nN+1) of (b_nN^k)
for n = -S.N_max:1:S.N_max
    Bessel_Evals(:, n+S.N_max+1) = ...
        besselj((n*S.N)+1, Bessel_Roots(:, n+S.N_max+1));
end

end