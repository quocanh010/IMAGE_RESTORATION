function output = imbilateral(input, varargin)
%% OUTPUT = IMBILATERAL(INPUT, TYPE, ...)
%
%  input/outputs:
%     INPUT:  input image
%     OUTPUT: filtered output
%
%  examples:
%     imbilateral(input, 'tau_space', radius, ...
%                 'tau_value', sensitivity, ...
%     imbilateral(input, 'tau_space', radius, ...
%                 'phi', phi, ...
%                 'central_weight', w)

import tools.*

options = makeoptions(varargin{:});

taus = getoptions(options, 'tau_space', [], 1);
if ~isfield(options, 'phi')
    tauv = getoptions(options, 'tau_value', [], 1);
    phi  = @(d) exp(-d.^2 / (2*tauv^2));
    w0   = exp(-1); % To avoid overweighting for the central pixel
else
    phi  = getoptions(options, 'phi', []);
    w0   = getoptions(options, 'central_weight', [], 1);
end

y = input;
x = w0 * input;
z = w0 * ones(size(input));
for k = -taus:taus
    for l = -taus:taus
        if (k == 0 && l == 0) || k^2 + l^2 >= taus^2
            continue;
        end
        y_shifted = imshift(y, k, l, varargin{:});

        % Exercise 6 (sheet #2)

        %----------------------------------------
        %----- FIXME: Complete this section -----
        %----------------------------------------
        distv = abs(y - y_shifted);
        %error('not implemented yet');
        %----------------------------------------
        
        w = phi(distv);

        %----------------------------------------
        %----- FIXME: Complete this section -----
        %----------------------------------------
        x = x + y_shifted .* w;
        %error('not implemented yet');
        %----------------------------------------


        z = z + w;
    end
end
output = x ./ z;
