function [output, info] = imconvolve(input, type, varargin)
%% OUTPUT = IMCONVOLVE(INPUT, TYPE, ...)
%
%  input/outputs:
%     INPUT:  input signal/image/video
%     OUTPUT: convolved signal/image/video
%     TYPE:   'gaussian', 'exponential', 'box', ...
%             'disk', 'laplacian', 'grad1', 'grad2', ...
%             'predefined'
%     INFO:   structure contining different information such as
%             the rate of noise reduction, the kernel function,
%             the fourier modulation...
%
%  examples:
%     imconvolve(input, 'box', 'tau', radius)
%     imconvolve(input, 'disk', 'tau', radius)
%     imconvolve(input, 'gaussian', 'tau', radius)
%     imconvolve(input, 'exponential', 'tau', radius)
%     imconvolve(input, 'exponential', 'rate', rate)
%     imconvolve(input, 'laplacian')
%     imconvolve(input, 'grad1')
%     imconvolve(input, 'grad2')
%     imconvolve(input, 'predefined', ...
%                'kernel', @kappa, 'tau', radius, ...
%                'normalize', boolean)
%     imconvolve(..., 'separable', boolean)
%     imconvolve(..., 'spectral',  boolean)
%     imconvolve(..., 'boundary', 'periodical'|'mirror'|'zeropadding'|'extension')


import tools.*

options = makeoptions(varargin{:});
spectral  = getoptions(options, 'spectral',  false);

% Define the radius 'tau' from the variance reduction factor 'rate'
if isfield(options, 'rate')
    if isfield(options, 'tau')
        warning('Field "rate" ignored: specify either "tau" or "rate"');
    else
        rate = getoptions(options, 'rate', [], 1);
        if rate < 1
            error('Rate of variance reduction must be bigger than 1');
        end
        switch type
            case {'gauss', 'gaussian'}
                tau = sqrt(rate / (4 * pi));
            case 'exponential'
                sumw  = @(tau) ((1+exp(-1./tau)) ./ (1-exp(-1./tau))).^2;
                sumw2 = @(tau) ((1+exp(-2./tau)) ./ (1-exp(-2./tau))).^2;
                tau   = fminsearch(@(tau) ...
                                   (sumw(abs(tau)).^2 ./ sumw2(abs(tau)) - rate).^2, ...
                                   sqrt(rate) / 4);
                tau   = abs(tau);
            case 'box'
                % Exercise 4 (sheet #5)

                %----------------------------------------
                %----- FIXME: Complete this section -----
                %----------------------------------------
                tau = (sqrt(rate) - 1)/2;
                %error('not implemented yet');
                %----------------------------------------

            case 'disk'
                tau = sqrt(rate / pi);
            otherwise
                error(sprintf('Rate option not supported for kernelt %s', type));
        end
        options.tau = tau;
    end
end

normalize = true;
switch type
    case {'gauss', 'gaussian'}
        separable = getoptions(options, 'separable', true);
        tau = getoptions(options, 'tau', 0, 1);
        hsize = ceil(5 * tau);
        kappa = @(i, j) exp(-(i.^2 + j.^2) / (2 * tau^2));
    case 'exponential'
        separable = getoptions(options, 'separable', false);
        tau = getoptions(options, 'tau', 0, 1);

        hsize = ceil(15 * tau);

        %----------------------------------------
        %----- DONE -----------------------------
        %----------------------------------------
        kappa = @(i,j) exp(-sqrt(i.^2 + j.^2) / tau);
        %----------------------------------------

    case 'box'
        separable = getoptions(options, 'separable', true);
        tau = getoptions(options, 'tau', 0, 1);
        hsize = ceil(tau);
        kappa = @(i, j) max(abs(i), abs(j)) <= tau;
    case 'disk'
        separable = getoptions(options, 'separable', false);
        tau = getoptions(options, 'tau', 0, 1);
        hsize = ceil(tau);


        %----------------------------------------
        %----- DONE -----------------------------
        %----------------------------------------
        kappa = @(i,j) i.^2 + j.^2 <= tau.^2;
        %----------------------------------------

    case 'laplacian'
        separable = getoptions(options, 'separable', false);
        normalize = false;
        hsize = 1;
        kappa = @(i, j) ...
                (i == 0 & j == 0) * -4 + ...
                (abs(i) + abs(j) == 1);
    case { 'grad1', 'grad1_forward' }
        separable = getoptions(options, 'separable', false);
        normalize = false;
        hsize = 1;
        kappa = @(i, j) ...
                (i == 0 & j == 0) * -1 + ...
                (i == 1 & j == 0);
    case { 'grad2', 'grad2_forward' }
        separable = getoptions(options, 'separable', false);
        normalize = false;
        hsize = 1;
        
        %----------------------------------------
        %----- DONE -----------------------------
        %----------------------------------------
        kappa = @(i, j) ...
                (i == 0 & j == 0) * -1 + ...
                (i == 0 & j == 1);
        %----------------------------------------

    case 'grad1_backward'
        separable = getoptions(options, 'separable', false);
        normalize = false;
        hsize = 1;
        kappa = @(i, j) ...
                (i == -1 & j == 0) * -1 + ...
                (i ==  0 & j == 0);
    case 'grad2_backward'
        separable = getoptions(options, 'separable', false);
        normalize = false;
        hsize = 1;
        kappa = @(i, j) ...
                (i == 0 & j == -1) * -1 + ...
                (i == 0 & j ==  0);
    case 'predefined'
        separable = getoptions(options, 'separable', false);
        normalize = getoptions(options, 'normalize', true);
        if ~spectral
            hsize = getoptions(options, 'tau', [], 1);
            kappa = getoptions(options, 'kernel', [], 1);
        else
            if isfield(options, 'kernel')
                kappa = getoptions(options, 'kernel', [], 1);
            else
                lambda = getoptions(options, 'modulation', [], 1);
            end
        end
    otherwise
        error(sprintf('convolution type "%s" not implemented', type));
end

normalization = 0;
if ~spectral
    output  = zeros(size(input));
    if ~separable
        % Exercise 5 (sheet #1)

        %----------------------------------------
        %----- DONE -----------------------------
        %----------------------------------------
        [n1, n2] = size(input);
        for k = -hsize:hsize
            for l = -hsize:hsize
                output = output +kappa(k,l) * imshift(input, -k, -l, varargin{:});
                normalization = normalization + kappa(k,l);
            end
        end
        %----------------------------------------

    else
        % Exercise 4 (sheet #2)

        %----------------------------------------
        %----- DONE -----------------------------
        %----------------------------------------
        normalization1 = 0;
        [n1, n2] = size(input);
        temp = output;
        for l = -hsize:hsize
            temp = temp +kappa(0,l) * imshift(input, 0, -l, varargin{:});
            normalization1 = normalization1 + kappa(0,l);
        end
        
        temp = temp / normalization1;
        for k = -hsize:hsize
            output = output +kappa(k,0) * imshift(temp, -k, 0, varargin{:});
            normalization = normalization + kappa(k,0);
        end
        %----------------------------------------


    end
else
    if isfield(options, 'boundary') && ~strcmp(options.boundary, 'periodical')
        error('Cannot use spectral modulation for non-periodical convolution');
    end

    % Define spectral modulations 'lambda' if not specified
    if ~exist('lambda', 'var')
        [n1, n2] = size(input);
        [u, v]   = fftgrid(n1, n2);

        % Exercise 3 (sheet #3)
        %----------------------------------------
        %----- DONE -----------------------------
        %----------------------------------------
        lambda = fft2(kappa(u,v));
        %----------------------------------------

    end

    % Exercise 3 (sheet #3)

    %----------------------------------------
    %----- DONE -----------------------------
    %----------------------------------------
    output = ifft2(lambda .* (fft2(input)));
    %----------------------------------------

    normalization = lambda(1, 1);
end

if normalize
    output = output / normalization;
end

% Return extra info if required
if nargout > 1
    % Desciption of the kernel
    info.type = type;
    if exist('kappa', 'var')
        info.kernel = kappa;
    end
    if exist('lambda', 'var')
        info.modulation = lambda;
    end
    if exist('tau', 'var')
        info.radius = tau;
    end
    if exist('hsize', 'var')
        info.hsize = hsize;
    end

    % Compute rate of variance reduction
    if ~spectral
        [i, j] = meshgrid(-hsize:hsize, -hsize:hsize);
        normalization = sum(sum(kappa(i, j)));

        % Exercise 4 (sheet #5)

        %----------------------------------------
        %----- FIXME: Complete this section -----
        %----------------------------------------
        %error('not implemented yet');
        rate = real(1./sum(sum (kappa(i,j).^2)));
        %----------------------------------------

    else
        rate = real(1 ./ mean(mean(abs(lambda).^2)));
    end
    if normalize
        rate = rate * normalization^2;
    end
    info.rate = rate;
end
