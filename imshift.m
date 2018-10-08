function output = imshift(input, k, l, varargin)
%% OUTPUT = IMSHIFT(INPUT, K, L, ...)
%
%  input/outputs:
%     INPUT:  input signal/image/video
%     OUTPUT: shifted signal/image/video
%     K, L:   displacement for first and second coordinate
%
%  examples:
%     imshift(input, 3, -10)
%     imshift(input, 5, 5, 'boundary', 'periodical'|'mirror'|'zeropadding'|'extension')

import tools.*

[n1, n2] = size(input);
options  = makeoptions(varargin{:});
boundary = getoptions(options, 'boundary', 'periodical');

% Fix bug: the shift was going to the wrong opposite sign
k = -k;
l = -l;

pk = max(0,-k);%100 
nk = max(0,+k);%0
pl = max(0,-l);%0
nl = max(0,+l);%50

switch boundary
    case 'periodical'
        xrange = mod((1:n1)-k - 1, n1)+1;% make an array mod from left to right then add 1
        yrange = mod((1:n2)-l - 1, n2)+1;
        output = input(xrange, yrange);
    case {'zeropadding', 'zero-padding', 'zero'}
        xrange = (pk+1):(n1-nk); %100 ->n1 // row
        yrange = (pl+1):(n2-nl); %1->n2-50 //collumn
        output = zeros(n1, n2);
        output(xrange+k, yrange+l) = input(xrange, yrange);
    case 'extension'
        % Exercise 5 (sheet #1)
%          xrange = (pk+1):(n1-nk);%101-> n1
%          yrange = (pl+1):(n2-nl);%1 -> n2 -50
%          output = zeros(n1, n2);
%          output(xrange+k, yrange+l) = input(xrange, yrange);
%          A = output(n1 - pk,yrange+nl);
%          C = repmat(A,100,1); 
%          output(((n1-99):n1),(nl+1:n2)) = C;
%          B = output(1:n1-pk,51);
%          D = repmat(B.',nl,1);
%          %F = output((1:n1-100),(1:50));
%          output((1:n1-100),(1:50)) = D.';
        xrange = mod((1:n1)-k - 1, n1)+1;% make an array mod from left to right then add 1
        yrange = mod((1:n2)-l - 1, n2)+1;
        xrange(n1-pk+1:n1) = n1;
        yrange(1:nl) = 1;
         yrange(n2-pl+1:n2) = n2;
         xrange(1:nk) = 1;
%         
        
%        xrange(1:nk) = pk+1;
%        xrange(n1-nk+1:n1) = n1-nk;
%yrange(1:nl) = pl+1;
       % yrange(n2-nl+1:n2) = n2-nl;
        %
       % yrange(1+pl:
       % output = input(xrange, yrange);
        %----------------------------------------
        %----- DONE -----------------------------
        %----------------------------------------

       
         output = input(xrange, yrange);

        %error('not implemented yet');

        %----------------------------------------

    case {'mirror', 'symmetrical'}
        % Exercise 5 (sheet #1)
        
        %----------------------------------------
        %----- DONE -----------------------------
        %----------------------------------------
        xrange = mod((1:n1)-k - 1, n1)+1;% make an array mod from left to right then add 1
        yrange = mod((1:n2)-l - 1, n2)+1;
        
        A = xrange(n1-pk:-1:n1-2*pk+1);
        xrange(n1-pk+1:n1) = A;
        xrange(1:nk) = xrange(2*nk-1:-1:nk);
        yrange(n2-pl:n2) = yrange(n2-pl:-1:n2-2*pl);
        %yrange(n2-pl:n2) = yrange(n2-pl:-1:n2-2*pl);
        yrange(1:nl) = yrange(2*nl:-1:nl+1);
        output = input(xrange, yrange);
        %----------------------------------------

    otherwise
        error(sprintf('boundary condition "%s" not implemented', boundary));
end
