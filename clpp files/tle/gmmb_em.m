%GMMB_EM     - EM estimated GMM parameters
%
% estS = gmmb_em(data)
% estS = gmmb_em(data, <params>...)
% [estS, stats] = gmmb_em(...)
%
% This version works with complex numbers too.
%
% data = N x D matrix
% params can be a list of 'name', value -pairs.
% stats is a matrix, row (cov fixes, loops, final log-likelihood)
%
% parameters (default value):
%
%   maxloops	maximum number of loops allowed (100)
%   thr		convergence threshold; log-likelihood change (1e-6)
%		set negative to use only maxloops condition
%   components	number of components in GMM (3)
%   verbose	print progress messages (false)
%
% Example:
%   estS = gmmb_em(data, 'components', 5, 'thr', 1e-8)
%
% References:
%   [1] Duda, R.O., Hart, P.E, Stork, D.G, Pattern Classification,
%   2nd ed., John Wiley & Sons, Inc., 2001.
%   [2] Bilmes, J.A., A Gentle Tutorial of the EM Algorithm and its
%    Application to Parameter Estimation for Gaussian Mixture and Hidden
%    Markov Models
%   International Computer Science Institute, 1998
%
% Author(s):
%    Joni Kamarainen <Joni.Kamarainen@lut.fi>
%    Pekka Paalanen <pekka.paalanen@lut.fi>
%
% Copyright:
%
%   Bayesian Classifier with Gaussian Mixture Model Pdf
%   functionality is Copyright (C) 2003 by Pekka Paalanen and
%   Joni-Kristian Kamarainen.
%
%   V_0_1 1.6  2004/02/19 06:45:01

function [estimate, varargout] = gmmb_em(data, varargin);

% default parameters
max_loops = 100;
em_threshold = 1e-6;
verbose = 0;
C = 3;
initS = [];

% parameter parsing
argc = nargin-1;
for argl = 1:2:(argc-1)
	switch lower(varargin{argl})
	  case 'maxloops'
	    max_loops = varargin{argl+1};
	  case 'thr'
	    em_threshold = varargin{argl+1};
	  case 'components'
	    C = varargin{argl+1};
	  case 'verbose'
	    verbose = varargin{argl+1};
	  case 'init'
	    initS = varargin{argl+1};
	  otherwise
%	    error(['unknown parameter ' varargin{argl}]);
	end
end



N = size(data,1);	% number of points
D = size(data,2);	% dimensions

if isempty(initS)
	% mus initialization (thanks V. Kyrki)
	if C>1
		mu = fcm(data, C, [2.0 100 1e-3 verbose]).';
		% fcm initialization has random nature, results will vary
	else
		mu = mean(data, 1).';
	end

	% covariances initialization
	nsigma = covfixer2(diag(diag(cov(data))));
	sigma = zeros(D,D,C);
	for c = 1:C
		sigma(:,:,c) = nsigma;
	end
else
	C = size(initS.mu, 2);
	mu = initS.mu;
	sigma = zeros(D,D,C);
	for c = 1:C
		sigma(:,:,c) = covfixer2(initS.sigma(:,:,c));
	end
end

% weights initialization
weight = ones(C,1) * (1/C);

% old values for stopping condition calculations
old_loglike = -realmax;

loops=1;
fixcount=0;

tulo = gmmcpdf(data, mu, sigma, weight);

while 1
	% one EM cycle
	pcompx = tulo ./ (sum(tulo,2)*ones(1,C));
	
	for c = 1:C
		% calculate new estimates
		psum = sum(pcompx(:,c));
		
		% weight
		weight(c) = 1/N*psum;
	
		% mean
		nmu = sum(data.*(pcompx(:,c)*ones(1,D)), 1).' ./ psum;
		mu(:,c) = nmu;
		
		% covariance
		moddata = (data - ones(N,1)*(nmu.')) .* (sqrt(pcompx(:,c))*ones(1,D));
		% sqrt(pcompx) is because it will be squared back
		nsigma = (moddata' * moddata) ./ psum;
		
		% covariance matrix goodness assurance
		sigma(:,:,c) = covfixer2(nsigma);
	end
	
	% finish test
	tulo = gmmcpdf(data, mu, sigma, weight);
	loglike = sum(log(sum(tulo, 2)));
	
	if verbose
		disp([ 'log-likelihood diff ' num2str(loglike-old_loglike)  ' on round ' num2str(loops) ]);
	end

	if abs(loglike - old_loglike) < em_threshold
		break;
	end
	
	if loops >= max_loops
		break;
	end

	loops = loops +1;
	old_loglike = loglike;
end



estimate = struct('mu', mu,...
		'sigma', sigma,...
		'weight', weight);

if(nargout>1)
	varargout(1) = {[fixcount loops loglike]};
end


% ------------------------------------------

function tulo = gmmcpdf(data, mu, sigma, weight);
N = size(data, 1);
C = size(weight,1);

pxcomp = zeros(N,C);
for c = 1:C
	pxcomp(:,c) = cmvnpdf(data, mu(:,c).', sigma(:,:,c));
end
tulo = pxcomp.*repmat(weight.', N,1);

