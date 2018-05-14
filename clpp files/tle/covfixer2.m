function nsigma = covfixer2(sigma);
%COVFIXER2   - force matrix to be a valid covariance matrix
%
% covmatrix = COVFIXER2(matrix)
% Matrix is forced (complex conjugate) symmetric,
% positive definite and its diagonal real valued.

%
% V_0_1
% covfixer2.m,v 1.7 2004/02/25 13:21:56 paalanen Exp
% Copyright 2003, Pekka Paalanen <pekka.paalanen@lut.fi>

% except isspd() function which is from The MathWorks Matlab mvnpdf.m.

D = size(sigma, 1);
covfixmat = ones(D) + 0.01*diag(ones(D,1));

nsigma = imagfixer(sigma);
while isspd(nsigma) == 0
	% covariance matrix is not positive definite
	% fix it
	if any(diag(nsigma) <= 0)
		% add 1% of max of diag to diag
		nsigma = nsigma + (covfixmat-1)*max(abs(diag(nsigma))+realmin);
	else
		% increase diagonal values by 1 percent
		nsigma = nsigma .* covfixmat;
	end
end


% ------------------

function [t,R] = isspd(Sigma)
%ISPDS Test if a matrix is positive definite symmetric
% T = ISPDS(SIGMA) returns a logical indicating whether the matrix SIGMA is
% square, symmetric, and positive definite, i.e., it is a valid full rank
% covariance matrix.
%
% [T,R] = ISPDS(SIGMA) returns the cholesky factor of SIGMA in R.  If SIGMA
% is not square symmetric, ISPDS returns [] in R.

%   Copyright 1993-2002 The MathWorks, Inc.
%   Revision: 1.2   Date: 2002/03/28 16:51:27

% Test for square, symmetric
[n,m] = size(Sigma);
if (n == m) & all(all(abs(Sigma - Sigma') < 10*eps*max(abs(diag(Sigma)))));
    % Test for positive definiteness
    [R,p] = chol(Sigma);
    if p == 0
        t = 1;
    else
        t = 0;
    end
else
    R = [];
    t = 0;
end

% ------------------

function nsigma = imagfixer(sigma);

% force symmetric
nsigma = sigma - (sigma - sigma')/2;
% purge imag
purge = imag(diag(nsigma));
nsigma = nsigma - diag(purge)*sqrt(-1);

if max(purge) > 1e-4
	WARNING('gmmbayes:covfixer2:imagfixer', 'Quite big imaginary components removed from the diagonal');
end
