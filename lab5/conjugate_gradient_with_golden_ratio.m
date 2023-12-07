% funkcja optymalizująca gradientami sprzężonymi z iteracyjnym złotym podziałem
function [x_opt, f_opt, n_iters] = conjugate_gradient_with_golden_ratio(ridgeFun, x0, max_iter, verbose)
    % ridgeFun: Function to minimize
    % x0: Initial guess for the parameters
    % max_iter: Maximum number of iterations
    
    tol = 1e-8;
    golder_ratio_tol = 1e-7;

    % Initialize parameters
    x = x0;
    [result, grad] = ridgeFun(x);
    p = -grad;
    
    if(verbose)
        fprintf('Iteration %d: f_opt = %f\n', 0, result);
    end

    % Stopping cryterium will be
    % when 25 times in a row there will be a small alpha
    max_time_small_alpha = 25;
    time_small_alpha = 0;
    small_alpha = 1e-7;

    % Stopping cryterium for times the function decreased only slightly
    max_time_decreased_slightly = 25;
    time_decreased_slightly = 0;
    decreased_slightly = 1e-4;

    % Perform conjugate gradient iterations
    for iter = 1:max_iter
        n_iters = iter;
        
        % Golden Ratio line search
        alpha = golden_ratio_search(@(a) ridgeFun(x + a * p), golder_ratio_tol);

        % Update parameters
        x_new = x + alpha * p;

        % Update gradient
        [result_new, grad_new] = ridgeFun(x_new);

        % Fast return in cast the numerical errors cause the function to
        % increase
        if (result_new > result)
            if verbose
                fprintf('numerical errors caused the function to increase in iteration %d\n', iter);
            end
            result_new = result;
            break;
        end

        % Compute beta for the next iteration
        beta = (grad_new' * grad_new) / (grad' * grad);

        % Update conjugate direction
        p = -grad_new + beta * p;

        % Update gradient for the next iteration
        grad = grad_new;
        x = x_new;

        % Display iteration information
        if verbose
            fprintf('Iteration %d: f_opt = %f\n', iter, result_new);
        end

        % stopping small slpha
        if (alpha < small_alpha)
            time_small_alpha = time_small_alpha + 1;
        else
            time_small_alpha = 0;
        end

        % stopping small improvement
        if (result - result_new < decreased_slightly)
            time_decreased_slightly = time_decreased_slightly + 1;
        else
            time_decreased_slightly = 0;
        end
        
        if (norm(grad) < tol) || (time_small_alpha > max_time_small_alpha) || (time_decreased_slightly > max_time_decreased_slightly)
            if verbose
                fprintf('Converged early at iteration %d\n', iter);
            end
            break;
        end

        result = result_new;
    end

    % Final result
    x_opt = x;
    f_opt = result_new;
end

function alpha = golden_ratio_search(fun, tol)
    % Golden Ratio line search
    alpha_low = 0;
    alpha_high = 1;

    golden_ratio = (sqrt(5) - 1) / 2;

    while abs(alpha_high - alpha_low) > tol
        alpha1 = alpha_low + (1 - golden_ratio) * (alpha_high - alpha_low);
        alpha2 = alpha_low + golden_ratio * (alpha_high - alpha_low);

        f1 = fun(alpha1);
        f2 = fun(alpha2);

        if f1 < f2
            alpha_high = alpha2;
        else
            alpha_low = alpha1;
        end
    end

    alpha = (alpha_low + alpha_high) / 2;
end
