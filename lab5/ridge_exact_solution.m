% funkcja znajdująca dokładnie rozwiązanie problemu Ridge analitycznie
function [x_opt, f_opt] = ridge_exact_solution(A, b, alpha)
    n = size(A, 2);
    x_opt = (A' * A + 2 * alpha * eye(n)) \ (A' * b);
    f_opt = (((A * x_opt - b)' * (A * x_opt - b)/2) + (alpha * x_opt' * x_opt)) / n;
end
