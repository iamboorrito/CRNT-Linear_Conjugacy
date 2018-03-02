m = 50;
n = 70;
ubound = 60;
ep = 0.8;

Y = randi(ubound-1, m, n);

A = randi(ubound-1, n, n);

for i=1:n
    A(i, i) = 0;
    A(i, i) = -sum(A(:, i));
end

M = Y*A;

fprintf('rows = %d;\ncols = %d;\nubound = %d;\neps = %f;\n',m, n, ubound, ep);

print_cplex_fmt('Y', Y)
print_cplex_fmt('M', M)
%print_cplex_fmt('A', A)
