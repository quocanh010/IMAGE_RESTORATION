function output = imnoise(input, type, varargin)
%% OUTPUT = IMNOISE(INPUT, TYPE, ...)
%
%  input/outputs:
%     INPUT:  input signal/image/video
%     OUTPUT: noisy signal/image/video
%     TYPE:   'gaussian', 'poisson', 'gamma', ...
%             'saltpepper', 'impulse'
%
%  examples:
%     imnoise(input, 'gaussian', 'sigma', value)
%     imnoise(input, 'poisson', 'Q', value)
%     imnoise(input, 'gamma', 'L', value)
%     imnoise(input, 'saltpepper', 'L', value, 'P', value)
%     imnoise(input, 'impulse', 'L', value, 'P', value)

import tools.*

options = makeoptions(varargin{:});
switch type
    case {'gauss', 'gaussian'}
        sigma = getoptions(options, 'sigma', 0, 1);
        output = input + sigma * randn(size(input));
    case {'shot', 'poisson', 'poissonian'}
        Q = getoptions(options, 'Q', 1);
        output = mypoissrnd(input / Q) * Q;
    case {'speckle', 'gamma'}
        L = getoptions(options, 'L', 1);
        output = mygamrnd(L, input / L);
    case {'rice', 'rician'}
        sigma = getoptions(options, 'sigma', 0, 1);
        theta = 2 * pi * rand(size(input));
        realp = input .* cos(theta) + sigma * randn(size(input))
        imagp = input .* sins(theta) + sigma * randn(size(input))
        output = sqrt(realp.^2 + imagp.^2);
    case {'saltpepper',  'impulse'}
        L = getoptions(options, 'L', 0, 1);
        P = getoptions(options, 'P', 0, 1);
        if sum(input(:) - round(input(:))) > 0 || ...
                min(input(:)) < 0 || max(input(:)) > L-1
            error(['input should only contain integer values in [0, L-1]']);
        end
        idx = rand(size(input)) < P;
        output = input;
        switch type
            case 'saltpepper'
                pattern = (L-1) * randi([0, 1], size(input));
            case 'impulse'
                pattern = randi([0, L-1], size(input));
        end
        output(idx) = pattern(idx);
    case 'ccd'
        Qe     = getoptions(options, 'Qe', 0, 1);
        D      = getoptions(options, 'D', 0, 1);
        t      = getoptions(options, 't', 0, 1);
        sigma  = getoptions(options, 'sigma', 0, 1);
        x      = input .* Qe .* t;
        l      = D .* t;
        z      = applynoise(x, 'poisson');
        n      = applynoise(l * ones(size(input)), 'poisson');
        w      = applynoise(zeros(size(input)), 'gaussian', 'sigma', sigma);
        output = z + n + w;
    otherwise
        error(sprintf('noise type "%s" not implemented', type));
end


function output = mypoissrnd(input)

if min(input(:)) < 0
    error(['input should only contain non-negative values']);
end
STEP = 500;
ll = input;
k = zeros(size(input));
p = ones(size(input));
nonstop = input < 50;
while sum(nonstop(:)) > 0
    k(nonstop) = k(nonstop) + 1;
    u = rand(size(input));
    p = p .* u;
    idx = nonstop & p < exp(1) & ll > 0;
    idx2 = idx & ll > STEP;
    idx3 = idx & ll <= STEP;
    p(idx2) = p(idx2) * exp(STEP);
    ll(idx2) = ll(idx2) - STEP;
    p(idx3) = p(idx3) .* exp(ll(idx3));
    ll(idx3) = -1;
    nonstop = nonstop & p > 1;
end
output = k - 1;
idx = input >= 50;
output(idx) = round(input(idx)) + sqrt(input(idx)) .* randn(sum(idx(:)), 1);

function output = mygamrnd(L, input)

if min(input(:)) < 0
    error(['input should only contain non-negative values']);
end
if L < 0
    error(['L must be a positive value']);
end

delta = L - floor(L);

nonstop = ones(size(input));
zeta = zeros(size(input));
eta = zeros(size(input));
while sum(nonstop(:)) > 0
    V3 = rand(size(input));
    V2 = rand(size(input));
    V1 = rand(size(input));
    v0 = exp(1) / (exp(1) + delta);
    idx1 = nonstop & V3 <= v0;
    idx2 = nonstop & V3 > v0;
    zeta(idx1) = V2(idx1).^(1 / delta);
    eta(idx1) = V1(idx1) .* zeta(idx1).^(delta-1);
    zeta(idx2) = 1 - log(V2(idx2));
    eta(idx2) = V1(idx2) .* exp(-zeta(idx2));
    nonstop = eta > zeta.^(delta-1) .* exp(-zeta);
end
output = zeros(size(input));
for k = 1:floor(L)
    output = output + log(rand(size(input)));
end
output = input .* (zeta - output);
