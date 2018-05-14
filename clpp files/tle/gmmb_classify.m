%GMMB_CLASSIFY Classify data using Bayesian or Mahalanobis distance classifier.
%
%     T = GMMB_CLASSIFY(S, data) Classifies D dimensional data (N points)
%     using Gaussian Mixture Model
%     Bayesian classifier in struct S into K classes.
%     S is a bayesS struct, see bayes_struct.txt.
%
%     See also GMMB_CREATE.
%
% Optional features:
%    T = GMMB_CLASSIFY(... , 'values')
%    Returns a posteriori class probabilities or Mahalanobis
%    distances instead of class labels.
%    T: N x K matrix
%
%    T = GMMB_CLASSIFY(... , 'mahalanobis')
%    Use Mahalanobis distances as classification rule
%    instead of Bayesian.
%
%
% References:
%   [1] Duda, R.O., Hart, P.E, Stork, D.G, Pattern Classification,
%   2nd ed., John Wiley & Sons, Inc., 2001.
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
%   V_0_1 1.6  2003/12/10 15:23:14
%

function [t] = gmmb_classify(bayesS, data_, varargin);

return_values = 0;
use_mahalanobis = 0;

% parameter parsing
argc = nargin-1;
for argl = 1:(argc-1)
	switch lower(varargin{argl})
	  case 'values'
	    return_values = 1;
	  case 'mahalanobis'
	    use_mahalanobis = 1;
	  otherwise
	end
end


N = size(data_,1);
K = size(bayesS,2);

P = [bayesS.apriories];

% data_ is N x D

if use_mahalanobis
	%disp('Using Mahalanobis distance.');
	% Mahalanobis distance classifier
	sqrmdist = zeros(N,K);
	for k = 1:K
		C = size(bayesS(k).mu, 2);
		sqrdist = zeros(N,C);
		for c = 1:C
			invs = inv(bayesS(k).sigma(:,:,c));
			mu = bayesS(k).mu(:,c).';
			sqrdist(:,c) = sum((data_*invs).*conj(data_),2) ...
				- data_*invs*mu' ...
				- (mu*invs*data_').' ...
				+ mu*invs*mu';
		end
		sqrmdist(:,k) = min(real(sqrdist), [], 2);
	end
	if return_values
		t = sqrmdist;
	else
		[a, b] = min(sqrmdist, [], 2);
		t = b';
	end
else
	% GMM Bayesian classifier
	
	% classify all points simultaneously
	pxomega = zeros(N,K);
	for k = 1:K
		pxomega(:,k) = gmmb_pdf(data_, bayesS(k).mu, bayesS(k).sigma, bayesS(k).weight );
	end
	tulo = pxomega.*repmat(P,N,1);
	PomegaKx = tulo ./ repmat(sum(tulo,2),1,K);

	if return_values
		t = PomegaKx;
	else
		[a, b] = max(PomegaKx , [], 2);
		t = b;
	end
end

