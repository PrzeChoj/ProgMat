%% test
x = randi([0 1500], n, 1);

ridgeFun(x)


%% yyy
n = 10;
m = 20;
alpha = 0.1;

[ridgeFun] = fun(n, m, alpha);


% Initial guess for beta
initialPoint = zeros(n, 1);

% Set options
options = optimset('Algorithm', 'quasi-newton', 'Display', 'iter', 'GradObj', 'on');
%options = optimset('Algorithm', 'quasi-newton', 'Display', 'iter');

% Use fminunc to find the argmin of ridgeFun with specified options
optimizedPoint = fminunc(ridgeFun, initialPoint, options);

disp('Optimized Coefficients:');
disp(optimizedPoint);
