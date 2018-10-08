function g = grad(x, varargin)
% Return the gradient field of an image

[n1, n2]   = size(x);
g          = zeros(n1, n2, 2);
g(:, :, 1) = imconvolve(x, 'grad1_forward', varargin{:});
g(:, :, 2) = imconvolve(x, 'grad2_forward', varargin{:});
