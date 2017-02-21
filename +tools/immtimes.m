function z = immtimes(x, y)
% Mulitply pointwise matrix field X with matrix filed Y 

[n1,  n2,  M, K] = size(x);
[n1b, n2b, Kb, N] = size(y);
if n1 ~= n1b || n2 ~= n2b || K ~= Kb
    error('Dimension mismatch');
end

z = zeros(n1, n2, M, N);
for i = 1:M
    for j = 1:N
        for k = 1:K
            z(:, :, i, j) = z(:, :, i, j) + ...
                x(:, :, i, k) .* y(:, :, k, j);
        end
    end
end
