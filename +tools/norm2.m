function n = norm2(x)
% Return the squared l2-norm of a vector field

n = sum(x.^2, 3);
