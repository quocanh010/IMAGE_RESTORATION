function options = makeoptions(varargin)

if mod(length(varargin), 2) == 1
    error('the number of options should be even');
end
options = struct();
for k = 1:2:length(varargin)
    options = setfield(options, varargin{floor(k/2)*2+1}, ...
                                varargin{floor(k/2)*2+2});
end
