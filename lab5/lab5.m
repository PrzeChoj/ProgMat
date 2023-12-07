% cd('/Users/adam/Desktop/programiki.nosync/ProgMat/lab5')
%% test
x0 = zeros(n, 1);

ridgeFun = fun(n, m, alpha);

%% test 2
[result, grad, hes] = ridgeFun(x0)

%% test 3
max_iter = 1000;
verbose = true;

[x_opt, f_opt] = conjugate_gradient_with_hessian(ridgeFun, x0, max_iter, verbose);
[x_opt_GR, f_opt_GR] = conjugate_gradient_with_golden_ratio(ridgeFun, x0, max_iter, verbose);
optimizedPoint = fminunc(ridgeFun, initialPoint, options);

x_opt_analitical = ridge_exact_solution(A, b, alpha);

sum(abs(x_opt - x_opt_analitical).^2)
sum(abs(optimizedPoint - x_opt_analitical).^2)
sum(abs(x_opt_GR - x_opt_analitical).^2)

[ridgeFun(x_opt) - ridgeFun(x_opt_analitical), ridgeFun(x_opt_GR) - ridgeFun(x_opt_analitical), ridgeFun(optimizedPoint) - ridgeFun(x_opt_analitical)]

%% etap 1
n = 3;
m = 5;
alpha = 0.1;

max_iter = 3;

[ridgeFun, A, b] = fun(n, m, alpha);


% Initial guess for beta
initialPoint = zeros(n, 1);

% Set options
options = optimset('Algorithm', 'quasi-newton', 'Display', 'iter', 'GradObj', 'on');

% Use fminunc to find the argmin of ridgeFun with specified options
optimizedPoint = fminunc(ridgeFun, initialPoint, options);

disp('Optimized Coefficients:');
disp(optimizedPoint);
