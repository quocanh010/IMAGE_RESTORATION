function h = plotimage(ima, range)

[n1, n2, K] = size(ima);
if sum(imag(ima(:))) == 0
    h = imagesc(ima);
    axis image;
    axis off;
    if exist('range', 'var')
        caxis(range);
    end
    if K == 1
        colormap(gray(1024));
    end
else
    import tools.quantile

    n_ima = abs(ima);
    a_ima = angle(ima);
    map = hsv(256);
    a = floor((a_ima + pi) / (2 * pi) * 255) + 1;
    qm = quantile(n_ima(:), 0.05);
    qM = quantile(n_ima(:), 0.95);
    w = (n_ima - qm) / (qM - qm);
    w(w > 1) = 1;
    w(w < 0) = 0;
    d1 = zeros(n1, n2);
    d2 = zeros(n1, n2);
    d3 = zeros(n1, n2);
    for k = 1:256
        idx = a == k;
        wi = w(idx);
        d1(idx) = (1 - wi) + wi * map(k, 1);
        d2(idx) = (1 - wi) + wi * map(k, 2);
        d3(idx) = (1 - wi) + wi * map(k, 3);
    end
    d = cat(3, d1, d2, d3);
    imagesc(d);
    %colormap(jet)
    axis image;
    axis off;
    if exist('range', 'var')
        caxis(range);
    end
end
if nargout < 1
    clear h;
end
