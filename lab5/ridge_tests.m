% funkcja testująca 100 razy
function [dokl, dokl_GR, dokl_MATLAB, n_iters, n_iters_GR, time_hess, time_GR, time_MATLAB] = ridge_tests(N, n, m, alpha, max_iter)
    n_iters = zeros(N, 1);
    n_iters_GR = zeros(N, 1);
    dokl = 0;
    dokl_GR = 0;
    dokl_MATLAB = 0;
    time_hess = 0;
    time_GR = 0;
    time_MATLAB = 0;

    for iter = 1:N
        [dokl_i, dokl_GR_i, dokl_MATLAB_i, n_iters_i, n_iters_GR_i, time_hess_i, time_GR_i, time_MATLAB_i] = single_test(n, m, alpha, max_iter);
        n_iters(iter) = n_iters_i;
        n_iters_GR(iter) = n_iters_GR_i;
        time_hess = time_hess + time_hess_i;
        time_GR = time_GR + time_GR_i;
        time_MATLAB = time_MATLAB + time_MATLAB_i;

        dokl = dokl + dokl_i;
        dokl_GR = dokl_GR + dokl_GR_i;
        dokl_MATLAB = dokl_MATLAB + dokl_MATLAB_i;
    end

    dokl = dokl / N;
    dokl_GR = dokl_GR / N;
    dokl_MATLAB = dokl_MATLAB / N;
end

function [dokl, dokl_GR, dokl_MATLAB, n_iters, n_iters_GR, time_hess, time_GR, time_MATLAB] = single_test(n, m, alpha, max_iter)
    tol = 1e-6;

    [ridgeFun, A, b] = fun(n, m, alpha);

    x0 = zeros(n, 1);
    verbose = false;
    options = optimset('Algorithm', 'quasi-newton', 'Display', 'off', 'GradObj', 'on');

    tic;
    [x_opt, ~, n_iters] = conjugate_gradient_with_hessian(ridgeFun, x0, max_iter, verbose);
    time_hess = toc;
    tic;
    [x_opt_GR, ~, n_iters_GR] = conjugate_gradient_with_golden_ratio(ridgeFun, x0, max_iter, verbose);
    time_GR = toc;
    tic;
    optimizedPoint = fminunc(ridgeFun, x0, options);
    time_MATLAB = toc;
    
    x_opt_analitical = ridge_exact_solution(A, b, alpha);

    dokl = sum(abs(x_opt - x_opt_analitical).^2);
    dokl_GR = sum(abs(x_opt_GR - x_opt_analitical).^2);
    dokl_MATLAB = sum(abs(optimizedPoint - x_opt_analitical).^2);
end