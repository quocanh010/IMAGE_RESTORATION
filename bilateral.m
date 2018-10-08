function output = bilateral(input,taus,sigma)


import tools.*

%% Parameters
%sigmaw_list = {10 20 30 60};
%taus        = 10;
tauv_factor = 3;

%% Core

    % Try to remove noise with bilateral filtering
    tauv = tauv_factor * sigma;
    output = imbilateral(abs(input), ...
                                 'tau_space', taus, ...
                                 'tau_value', tauv, ...
                                 'boundary', 'mirror');

end