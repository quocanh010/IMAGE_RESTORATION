function data = loaddata()

f         = fopen('./data/Abdomen.dat', 'r');
data      = double(fread(f, 'int16'));
dimx      = 128;
dimy      = 128;
dymt      = 200;
data      = reshape(data, dimx, dimy, dymt, 2);
data      = data(:,:,:,1) .* exp(i * data(:,:,:,2)/4096 * 2 * pi);
data      = permute(data, [2 1 3]);
