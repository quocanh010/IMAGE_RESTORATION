function q = quantile(x, p)

n = length(x);
x = sort(x, 'ascend');
q = x(ceil((n-1) * p) + 1);
