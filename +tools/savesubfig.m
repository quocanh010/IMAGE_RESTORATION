function savesubfig(h, dirname)
%% SAVESUBFIG(H, DIRNAME)
%
%  save the figure h in directory DIRNAME and all its subplots as
%    DIRNAME/main.fig
%    DIRNAME/subplot1.fig
%    DIRNAME/subplot2.fig
%    ...
%
%  for subaxes containing image data, SAVEFIGSUB will also create
%    DIRNAME/subimg1.png
%    DIRNAME/subimg2.png
%    ...
%
%  Copyright 2017 Charles Deledalle

[status, msg,~] = mkdir(dirname);
if ~status
    error(sprintf('Cannot create directory "%s": %s', dirname, msg));
end
savefig(h, [ dirname '/' 'main.fig' ]);
haxes = findobj(h, 'type', 'axes');
K = length(haxes);
for k = 1:K
    haxe = haxes(K-k+1);
    [img, flag] = getimage(haxe);
    switch flag
        case 0 % Not an image
        case 1 % Indexed image
            warning('case 1 not implemented yet')
        case 2 % Intensity image in standard range
            imwrite(img, [ dirname '/' sprintf('subimg%d.png', k) ]);
        case 3 % Intensity image not in standard range
            switch get(haxe, 'climmode')
                case 'manual'
                    range = get(haxe, 'clim');
                case 'auto'
                    range = [min(img(:)) max(img(:))];
            end
            img = (img - range(1)) / (range(2) - range(1));
            imwrite(img, [ dirname '/' sprintf('subimg%d.png', k) ]);
        case 4 % RGB Image
            imwrite(img, [ dirname '/' sprintf('subimg%d.png', k) ]);
        case 5 % Binary image
            warning('case 5 not implemented yet')
    end
    hsubfig = figure('visible', 'off');
    set(hsubfig,'CreateFcn','set(gcf,''Visible'',''on'')')
    haxe_new = copyobj(haxe, hsubfig);
    set(haxe_new, 'Position', get(0, 'DefaultAxesPosition'));
    colormap(hsubfig, colormap(h));
    filename = [ dirname '/' sprintf('subplot%d.fig', k) ];
    savefig(hsubfig, filename);
    delete(hsubfig);
end
