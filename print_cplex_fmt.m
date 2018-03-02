function [ ] = print_cplex_fmt(name, A )

[m, n] = size(A);

fprintf('%s = #[\n', name);
for i = 1:m
    fprintf('%d: [', i);
    
    for j = 1:n-1
        fprintf('%d, ', A(i, j));
    end
    
    fprintf('%d]', A(i, n));
    
    if i < m
        fprintf(',\n');
    end
end
fprintf('\n]#;\n');

end

