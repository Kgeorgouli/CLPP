%GMMB_GEM    - Greedy EM estimated GMM parameters
% Produces a bayesS struct without 'apriories'
% This is just a wrapper for the Vlassis Greedy EM algorithm implementation.
%
% estimate = GMMB_GEM(data[, parameters])
% Parameters (default):
%   verbose	print some progress numbers (false)
%   animate	plot data and ellipses during algorithm evaluation (false)
%   Cmax	the maximum number of GMM components
%   ncand	number of candidate locations for each new component (10)
% At least Cmax should be set explicitly.
% example:
%    estS = gmmb_gem(data, 'Cmax', 10, 'animate', true);
%
% References:
%   [1] Vlassis, N., Likas, A., A Greedy EM Algorithm for Gaussian Mixture
%   Learning, Neural Processing Letters 15, Kluwer Academic Publishers, 2002.
%   http://carol.wins.uva.nl/~vlassis/research/learning/index_en.html
%
% Author(s):
%    Pekka Paalanen <pekka.paalanen@lut.fi>
%
% Copyright:
%
%   Bayesian Classifier with Gaussian Mixture Model Pdf
%   functionality is Copyright (C) 2003 by Pekka Paalanen and
%   Joni-Kristian Kamarainen.
%
%   V_0_1 1.4  2004/02/19 12:59:25
%

function [estimate, varargout] = gmmb_gem(data, varargin);

[N, D] = size(data);	% number of points (n), dimensions (d)

% defaults
C = ceil(min(50, N/(D*D)/3));	% max n of components (k_max)
verbose = 0;
ncand = 10;
animation = 0;

% parameter parsing
argc = nargin-1;
for argl = 1:2:(argc-1)
	switch lower(varargin{argl})
	  case 'cmax'
	    C = varargin{argl+1};
	  case 'verbose'
	    verbose = varargin{argl+1};
	  case 'ncand'
	    ncand = varargin{argl+1};
	  case 'animate'
	    animation = varargin{argl+1};
	  otherwise
	end
end

[W, M, R, Tlogl] = gmmbvl_em(data, [], C, ncand, animation, verbose);

Cfinal = size(R,1);
sigma = zeros(D, D, Cfinal);
for c = 1:Cfinal
	Rk = reshape(R(c,:),D,D);
	sigma(:,:,c) = Rk' * Rk;
end

estimate = struct('mu', M.',...
	'sigma', sigma,...
	'weight', W);

if(nargout>1)
	varargout(1) = {Tlogl};
end
