function output = imdiffuse(input, type, varargin)
%% OUTPUT = IMDIFFUSE(INPUT, TYPE, ...)
%
%  input/outputs:
%     INPUT:  input signal/image/video
%     OUTPUT: diffused signal/image/video
%     TYPE:   'heat': Heat equation,
%             'oad':  Original anisotropic diffusion,
%             'rad':  Regularized anisotropic diffusion,
%             'tad':  Truly anisotropic diffusion,
%             'predefined': Solve dx/dt = A(x, x)
%                           with A(x, z) = A(x) z
%
%  examples:
%     imdiffuse(input, type, 'time', t)
%     imdiffuse(input, type, 'nsteps', m)
%     imdiffuse(input, type, ..., 'alpha', alpha = 1);
%     imdiffuse(input, type, ..., 'gamma', gamma = 0.99)
%     imdiffuse(input, type, ..., 'tau',   tau = 0.249)
%     imdiffuse(input, 'predefined, 'A', @(x, z) A(x, z))

import tools.*

options = makeoptions(varargin{:});

% Setting the step parameters
[n1, n2]   = size(input);
n          = n1 * n2;
dx2        = 1 / n;
alpha      = getoptions(options, 'alpha', 1);
gamma      = getoptions(options, 'gamma', 0.99);
if isfield(options, 'time')
    t      = getoptions(options, 'time', [], 1);
    if isfield(options, 'nsteps')
        m  = getoptions(options, 'nsteps', [], 1);
    else
        m  = ceil(4 * alpha * t / dx2 / gamma);
    end
    dt     = t / m;
else
    m      = getoptions(options, 'nsteps', [], 1);
    dt     = gamma * dx2 / (4 * alpha);
end
tau        = getoptions(options, 'tau', alpha * dt / dx2);
boundary   = getoptions(options, 'boundary', 'periodical');
scheme     = getoptions(options, 'scheme', 'explicit');

% Setting the g function for anisotropic diffusion
if ~strcmp(type, 'heat') && ~strcmp(type, 'predefined')
    if ~isfield(options, 'g')
        beta   = getoptions(options, ...
                            'beta', quantile(vect(norm2(grad(input))), 0.04));
        g      = @(d) alpha * exp(-d / beta);
    else
        g      = getoptions(options, 'g', []);
    end
end

% Setting the type of PDE as: dx/dt = A(x, x)
% where A(x, z) = A(x) z
switch type
    case 'heat'
        A = @(x, z) imconvolve(z, 'laplacian', varargin{:});
    case 'oad'
        A = @(x, z) A_oad(x, z, g, varargin{:});
    case 'rad'
        sigma = getoptions(options, 'sigma', 1);
        A = @(x, z) A_rad(x, z, g, sigma, varargin{:});
    case 'tad'
        sigma = getoptions(options, 'sigma', 1);
        rho   = getoptions(options, 'rho', 0.5);
        A = @(x, z) A_tad(x, z, g, sigma, rho, varargin{:});
    case 'predefined'
        A = getoptions(options, 'A');
    otherwise
        error(sprintf('type %s not implemented', type));
end

% Define laplacian kernel if required
if strcmp(type, 'heat') && strcmp(boundary, 'periodical') && ~strcmp(scheme, 'continuous')
    L = zeros(n1, n2);
    L(1, 1)  = -4;
    L(1, 2)  = 1;
    L(2, 1)  = 1;
    L(1, n2) = 1;
    L(n1, 1) = 1;
    L = fft2(L);
end

% Core
x = input;
switch scheme
    case 'solution'
        if strcmp(type, 'heat') && strcmp(boundary, 'periodical')
            [u, v]   = fftgrid(n1, n2);
            h        = sqrt(2 * tau * m);
            K        = exp(-(u.^2 + v.^2) / (2 * h^2)) / (2 * pi * h^2);
            K        = fft2(K);

            % Exercise 3 (sheet #4)

            %----------------------------------------
            %----- FIXME: Complete this section -----
            %----------------------------------------
            x = real(ifft2(fft2(x) .* K));
            %error('not implemented yet');
            %----------------------------------------

        else
            warning('solution not known in closed form: returned NaN');
            x = nan * ones(n1, n2);
        end
    case 'explicit'
        if strcmp(type, 'heat') && strcmp(boundary, 'periodical')
            % Exercise 3 (sheet #4)

            %----------------------------------------
            %----- FIXME: Complete this section -----
            %----------------------------------------
            K = (1+tau .*L).^m;
            %x = real(ifft2(fft2(x) .* Kee));
            %error('not implemented yet');
            %----------------------------------------

            x = real(ifft2(fft2(x) .* K));
        else
            for k = 1:m
                x = x + tau * A(x, x);
            end
        end
    case 'implicit'
        if strcmp(type, 'heat') && strcmp(boundary, 'periodical')
            K = 1 ./ (1 - tau .* L).^m;
            x = real(ifft2(fft2(x) .* K));
        else
            if strcmp(boundary, 'periodical')
                for k = 1:m
                    % Exercise 4 (sheet #4)

                    %----------------------------------------
                    %----- FIXME: Complete this section -----
                    %----------------------------------------
                    x = imcgs (@(z) z - tau * A(x,z));
                    %error('not implemented yet');
                    %----------------------------------------

                end
            else
                for k = 1:m
                    x = imgmres(@(z) z - tau * A(x, z), x);
                end
            end
        end
    otherwise
        error(sprintf('scheme %s not implemented', scheme));
end
output = x;



%%% Original Anisotropic Diffusion
function A = A_oad(x, z, g, varargin)

import tools.*

a  = g(norm2(grad(x, varargin{:})));
A  = div(bsxfun(@times, a, grad(z, varargin{:})), varargin{:});



%%% Regularized Anisotropic Diffusion
function A = A_rad(x, z, g, sigma, varargin)

import tools.*

x_conv = imconvolve(x, 'gaussian', 'tau', sigma, varargin{:});
% Exercise 4 (sheet #4)

%----------------------------------------
%----- FIXME: Complete this section -----
%----------------------------------------
a  = g(norm2(grad(x_conv, varargin{:})));
%error('not implemented yet');

%----------------------------------------

A = div(bsxfun(@times, a, grad(z, varargin{:})), varargin{:});



%%% Truly Anisotropic Diffusion
function A = A_tad(x, z, g, sigma, rho, varargin)

import tools.*

[n1, n2] = size(x);
x_conv   = imconvolve(x, 'gaussian', 'tau', sigma, varargin{:});
g_conv   = grad(x_conv);
M_conv   = zeros(n1, n2, 2, 2); % Tensor field capturing the average
                                % direction of the gradient
for k = 1:2
    for l = 1:2
        % Exercise 4 (sheet #4)

        %----------------------------------------
        %----- FIXME: Complete this section -----
        %----------------------------------------
        M_conv(:, :, k, l) = (g_conv(:,:,k) ) .* (g_conv(:,:,l));
        %error('not implemented yet');
        %----------------------------------------

        M_conv(:, :, k, l) = imconvolve(M_conv(:, :, k, l), ...
                                        'gaussian', 'tau', rho, ...
                                        varargin{:});
    end
end
T = imspecfunc(M_conv, g);      % Apply g to the eigenvalues of each tensor
A = div(immtimes(T, grad(z, varargin{:})), varargin{:});
