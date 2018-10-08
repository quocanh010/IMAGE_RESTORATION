function output = imosf(input, type, varargin)
%% OUTPUT = IMOSF(INPUT, TYPE, ...)
%
%  input/outputs:
%     INPUT:  input image
%     OUTPUT: filtered output
%     TYPE:   'gaussian', 'exponential', 'box', ...
%             'disk', 'laplacian', 'grad1', 'grad2', ...
%             'predefined'
%
%  examples:
%     imosf(input, 'median',   'tau', radius)
%     imosf(input, 'erode',    'tau', radius)
%     imosf(input, 'dilate',   'tau', radius)
%     imosf(input, 'open',     'tau', radius)
%     imosf(input, 'close',    'tau', radius)
%     imosf(input, 'trimmed',  'tau', radius, 'keep', k)
%     imosf(input, 'extremal', 'tau', radius)

import tools.*

options = makeoptions(varargin{:});
hsize = getoptions(options, 'tau', [], 1);
n3  = (2 * hsize + 1)^2;
switch type
    case 'median'
        filter = @(x) x(:, :, ceil(n3 / 2));
    case 'erode'
        filter = @(x) x(:, :, 1);
    case 'dilate'
        % Exercise 4 (sheet #2)

        %----------------------------------------
        %----- FIXME: Complete this section -----
        %----------------------------------------
        filter = @(x) x(:, :, n3);
        %error('not implemented yet');
        %----------------------------------------

    case 'open'
        output = imosf(input,  'erode', varargin{:});
        output = imosf(output, 'dilate', varargin{:});
        return
    case 'close'
        % Exercise 4 (sheet #2)

        %----------------------------------------
        %----- FIXME: Complete this section -----
        %----------------------------------------
        output = imosf(input,  'dilate', varargin{:});
        output = imosf(output, 'erode', varargin{:});
        return
        
        %----------------------------------------

    case 'trimmed'
        k = getoptions(options, 'keep', ceil(n3/2));
        filter = @(x) mean(trimmed(x, k), 3);
    case 'extremal'
        filter = @(x) extremal(x);
end

%% OSF filtering (only based on ranking)
[n1, n2] = size(input);
stack = zeros(n1, n2, n3);
c = 1;
for k = -hsize:hsize
    for l = -hsize:hsize
        stack(:, :, c) = imshift(input, k, l, varargin{:});
        c = c + 1;
    end
end
stack  = sort(stack, 3, 'ascend');
output = filter(stack);



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function x = trimmed(x, k)
% Keep the k central elements of x in the third dimension

n = size(x, 3);
if k <= 0
    error(sprintf('cannot keep %d elements', k))
end
if k > n
    k = n;
end
r = n - k;
x = x(:, :, (1+ceil(r/2)):(n-(r-ceil(r/2))));



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function output = extremal(stack)

[n1, n2, n3] = size(stack);
omin  = stack(:, :, 1);
omax  = stack(:, :, n3);
omean = mean(stack, 3);
output = zeros(n1, n2);

% Exercise 4 (sheet #2)

%----------------------------------------
%----- FIXME: Complete this section -----
%----------------------------------------
disi = abs( omean- omin);
disa = abs( omean -omax);
omin(disi > disa) = 0;
omax(disa > disi) = 0;
output = omin + omax;

    
%error('not implemented yet');
%----------------------------------------

