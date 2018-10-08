function output = median_filter(input, radius)
import tools.*

sigmaw = 20;
[n1, n2] = size(input);
output = imosf(abs(input), 'median', 'tau', radius);
end