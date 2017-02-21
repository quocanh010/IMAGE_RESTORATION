close all;
clear all;

import tools.*

% Load data
data = loaddata();
[n1, n2, T] = size(data);

% Display data
h = fancyfigure;
for k = 1:T
    subplot(131);
    plotimage(abs(data(:, :, k)), [0 255]);
    title('Magnitude image');
    subplot(132);
    plotimage(angle(data(:, :, k)));
    title('Phase image');
    subplot(133);
    plotimage(data(:, :, k));
    title('Complex image');
    drawnow
    pause(0.01);
end
