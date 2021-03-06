function d = KG_clppL2( a, b,isdiag )
% L2_DISTANCE - computes Euclidean distance matrix
%
% E = L2_distance(A,B)
%
%    A - (DxM) matrix 
%    B - (DxN) matrix
% 
% Returns:
%    E - (MxN) Euclidean distances between vectors in A and B
%
%
% Description : 
%    This fully vectorized (VERY FAST!) m-file computes the 
%    Euclidean distance between two vectors by:
%
%                 ||A-B|| = sqrt ( ||A||^2 + ||B||^2 - 2*A.B )
%
% Example : 
%    A = rand(400,100); B = rand(400,200);
%    d = distance(A,B);

% Author   : Roland Bunschoten
%            University of Amsterdam
%            Intelligent Autonomous Systems (IAS) group
%            Kruislaan 403  1098 SJ Amsterdam
%            tel.(+31)20-5257524
%            bunschot@wins.uva.nl
% Last Rev : Wed Oct 20 08:58:08 MET DST 1999
% Tested   : PC Matlab v5.2 and Solaris Matlab v5.3

% Copyright notice: You are free to modify, extend and distribute 
%    this code granted that the author of the original code is 
%    mentioned as the original author of the code.

% Fixed by JBT (3/18/00) to work for 1-dimensional vectors
% and to warn for imaginary numbers.  Also ensures that 
% output is all real, and allows the option of forcing diagonals to
% be zero.  
%
%

% This file is part of the Matlab Toolbox for Dimensionality Reduction v0.4b.
% The toolbox can be obtained from http://www.cs.unimaas.nl/l.vandermaaten
% You are free to use, change, or redistribute this code in any way you
% want for non-commercial purposes. However, it is appreciated if you 
% maintain the name of the original author.
%
% (C) Laurens van der Maaten
% Maastricht University, 2007

    if(~exist('isdiag','var'))
        isdiag=0;
    end
    if nargin < 2
       error('Not enough input arguments');
    end
    if size(a, 1) ~= size(b, 1)
        error('A and B should be of same dimensionality');
    end
    if ~isreal(a) || ~isreal(b)
        warning('Computing distance table using imaginary inputs. Results may be off.'); 
    end

    % Padd zeros if necessray
    if size(a, 1) == 1
        a = [a; zeros(1, size(a, 2))]; 
        b = [b; zeros(1, size(b, 2))]; 
    end

    if(size(a,2)<8000)
        % Compute distance table
        aa = sum(a .* a);
        bb = sum(b .* b);
        ab = a' * b;
        d = sqrt(repmat(aa', [1 size(bb, 2)]) + repmat(bb, [size(aa, 2) 1]) - 2 * ab);
    else
        warning('matrices too large???');
        d=zeros(size(a,2),size(b,2));
        chunkSize=6000;
        blockNums = ceil(size(b,2)/chunkSize);
        aa = sum(a .* a);
        for i=1:blockNums
            maxB = i*chunkSize;
            if(maxB>size(b,2))
                maxB=size(b,2);
            end
            bbb = b(:,(i-1)*chunkSize+1:maxB);
            bb = sum(bbb .* bbb);
            ab = a' * bbb;
            d(:,(i-1)*chunkSize+1:maxB) = sqrt(repmat(aa', [1 size(bb, 2)]) + repmat(bb, [size(aa, 2) 1]) - 2 * ab);
        end
    end

    % Make sure result is real
    d = real(d);

    if any(any(d<0))
        d(d<0) = 0;
    end

    if(isdiag==1)
        %enforce 0 on diagonal
        d = d.*(1-eye(size(d)));
    end

end

