% funkcja optymalizująca gradientami sprzężonymi
function [x_opt, f_opt, n_iters] = conjugate_gradient_with_hessian(ridgeFun, x0, max_iter, verbose)
    % ridgeFun: Function to minimize
    % x0: Initial guess for the parameters
    % max_iter: Maximum number of iterations

    tol = 1e-8;
    learningRate = 1;

    % Initialize parameters
    x = x0;
    [result, grad, hes] = ridgeFun(x);
    p = -grad;

    if(verbose)
        fprintf('Iteration %d: f_opt = %f\n', 0, result);
    end

    % Perform conjugate gradient iterations
    for iter = 1:max_iter
        n_iters = iter;
        
        % Update parameters using the Newton method
        x_new = x - learningRate * (hes \ grad);

        % Update gradient and hessian
        [result_new, grad_new, hes_new] = ridgeFun(x_new);

        % Fast return in cast the numerical errors cause the function to
        % increase
        if (result_new > result)
            result_new = result;
            break;
        end
        
        % Compute beta for the next iteration
        beta = (grad_new' * grad_new) / (grad' * grad);

        % Update conjugate direction
        p = -grad_new + beta * p;

        % Update gradient for the next iteration
        grad = grad_new;
        hes = hes_new;
        x = x_new;

        % Display iteration information
        if(verbose)
            fprintf('Iteration %d: f_opt = %f\n', iter, result_new);
        end

        if norm(grad) < tol
            if(verbose)
                fprintf('Converged early at iteration %d\n', iter);
            end
            break;
        end
    end

    % Final result
    x_opt = x;
    f_opt = result_new;
end
