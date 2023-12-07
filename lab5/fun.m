% Zwraca funkcję do robienia Ridge oraz jej gradient oraz Hesian
function [ridgeFun] = fun(n, m, alpha)
    [A, b] = drawData(n, m);

    ridgeFun = @innerFunction;

    function [result, grad, hes] = innerFunction(x)
        result = ((A * x - b)' * (A * x - b)/2) + (alpha * x' * x);
        grad = A' * (A * x - b) + 2 * alpha * x;
        hes = A' * A + 2 * alpha * eye(n);
    end
end


% Funkcja drawData generuje dane do testów.
function [A, b] = drawData(n, m)
    A = randi([0 1500], m, n);
    b = randi([0 1500], m, 1);
end

