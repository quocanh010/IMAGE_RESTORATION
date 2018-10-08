clear all
close all

import tools.*
% Assume that the noise is constant over time.
%% Parameters

homomorphic  = false;
anscombe     = false;

%% Core

data = loaddata();
[n1, n2, K] = size(data);


%% Exercise 3 (sheet #1)
if homomorphic
     %y_gaussian = log(y_gaussian);
    
     %y_poisson = log(y_poisson);
    %----------------------------------------
    %----- FIXME: Complete this section -----
    %----------------------------------------
    y_gamma = log(y_gamma);
    %----------------------------------------

end
if anscombe
   
    %----------------------------------------
    %----- FIXME: Complete this section -----
    %----------------------------------------
    y_poisson = 2*sqrt(y_poisson + 3/8);
    %----------------------------------------

end

%% Statistics
my_mean = mean(abs(data), 3);
sy_var = std(abs(data), [], 3);

%% Display
h = fancyfigure;
subplot(2,1,1);
plot(my_mean(:), sy_var(:), 'x');
xlim([0, 256])
xlabel('$mean$');
ylabel('$var$');
axis square
subplot(2,1,1);


