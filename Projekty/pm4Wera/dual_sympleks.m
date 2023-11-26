function [ZPx, ZDy, exitflag] = dual_sympleks(c, A, b, g)
c
A
b
g
exitflag =1;
[n_rows, n_cols] = size(A');
A_dual = [A' eye(n_rows) -eye(n_rows)];
c_dual = [b' g' zeros(1, n_rows)];
b_dual = c;
base_indices = zeros(n_rows, 1);
sign_changed = zeros(n_rows, 1);
A_dual
b_dual
c_dual
for i=1:n_rows
    if(b_dual(i)<0)
        sign_changed(i) = 1;
        base_indices(i) = n_cols + n_rows + i;
        A_dual(i, :) = -A_dual(i, :);
        b_dual(i, :) = -b_dual(i, :);
    else
        base_indices(i) = n_cols + i;
        
    end
end
sign_changed
A_dual
b_dual
c_dual
lb = zeros(size(b_dual))';
[ZDy, ~, exitflag, A_dual, ~, indices] = sympleks(-c_dual, A_dual, b_dual, base_indices);
linprog(c_dual, [], [], A_dual, b_dual, lb)
ZDy
A_solution = A_dual(:, base_indices);
for i=1:n_rows
    if(sign_changed(i)==1)
        A_solution(:, i) = -A_solution(:, i);
        %A_dual(:, base_indices(i)) = -A_dual(:, base_indices(i));
    end
end
ZPx = c_dual(indices)* A_solution;
end