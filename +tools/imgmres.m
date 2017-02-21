function [c, flag] = imgmres(A, b, varargin)
% Solve A b = c with A square
% with Generalized Minimum Residual Method
% Unlike gmres, b, c does not need to be vectors

mat       = @(x) reshape(x, size(b));
vect      = @(x) x(:);
[c, flag] = gmres( @(x) vect(A(mat(x))), vect(b), varargin{:});
c         = mat(c);
