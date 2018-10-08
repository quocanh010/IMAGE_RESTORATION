function d = div(x, varargin)
% Return the divergence of a vector field

d = imconvolve(x(:, :, 1), 'grad1_backward', varargin{:}) + ...
    imconvolve(x(:, :, 2), 'grad2_backward', varargin{:});
