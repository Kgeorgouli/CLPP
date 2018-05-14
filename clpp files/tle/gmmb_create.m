%GMMB_CREATE - Construct new Bayesian classifier with Gaussian mixture model pdf
%
%     S = GMMB_CREATE(data, class, method [, parameters]) Generates a Bayesian
%     classifier for two classes both having GMM distributed pdf
%     with estimated mean values, variances and
%     apriories. Classifier is returned in bayesS struct S.
%     method can be 'EM', 'FJ' or 'GEM'.
%     EM and FJ can work with complex numbers.
%
%     See also GMMB_CLASSIFY.
%
%     Parameters are delegated directly to underlaying GMM estimation
%     function (gmmb_em, gmmb_fj, gmmb_gem). See also them.
%
% Examples:
%   bayesS = gmmb_create(data, class, 'EM', 'components', 5, 'thr', 1e-8);
%   bayesS = gmmb_create(data, class, 'FJ', 'Cmax', 50, 'thr', 1e-9);
%   bayesS = gmmb_create(data, class, 'GEM', 'Cmax', 10, 'verbose', true);
%
%   The bayesS struct is documented in bayes_struct.txt.
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
%   Bayesian Classifier with Gaussian Mixture Model Pdf is
%   Copyright (C) 2003 by Pekka Paalanen and Joni-Kristian
%   Kamarainen.
%
%   V_0_1 1.4  2003/08/29 10:48:17
%

function bayesS = gmmb_create(data, cl, method, varargin);

if strcmp(method, 'SCC_EM')
	% SCC_EM is an "undocumented" method, SCC initialized EM.
	% It requires the drfwnewscc.m from the facedetect module.
	%sccS = drfwnewscc(data, cl, 'NC', 14);
	configure;
	sccS = load(sccClassifierSaveFile, 'classifierS');
	sccS = sccS.classifierS;
end

K = max(cl);

mu ={};
sigma = {};
weight = {};
prior = {};

for k = 1:K
	cvals = data(cl == k, :);
	N = size(cvals,1);	% points
	if N<2
		error('MATLAB:badopt', ['not enough data points in training set #' num2str(k)]);
	end
	
	switch method
	case 'EM'
		estim = gmmb_em(cvals, varargin{:});
	case 'SCC_EM'
		init_sigma = shiftdim(sccS.sampleCovs,1);
		initS = struct(...
			'sigma', init_sigma(:,:,(sccS.classInds==k)), ...
			'mu', sccS.classMeans((sccS.classInds==k), :).' );
		estim = gmmb_em(cvals, 'init', initS, varargin{:});
	case 'FJ'
		estim = gmmb_fj(cvals, varargin{:});
	case 'GEM'
		estim = gmmb_gem(cvals, varargin{:});
	otherwise
		error('Unknown method');
	end

	mu{k} = estim.mu;
	sigma{k} = estim.sigma;
	weight{k} = estim.weight;
	prior{k} = N/size(data,1);
end

bayesS = struct('mu', mu,...
		'sigma', sigma,...
		'apriories', prior,...
		'weight', weight);
