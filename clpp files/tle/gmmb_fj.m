%GMMB_FJ     - Figueiredo-Jain estimated GMM parameters
% Produces a bayesS struct without 'apriories' or NaN in case of failure.
%
% Works with complex numbers directly.
%
% estimate = GMMB_FJ(data[, parameters])
% Parameters (default):
%   Cmax	the maximum number of GMM components
%		set to -1 to use all data points as component means
%   Cmin	the minimum number of GMM components tried (1)
%   verbose	print some progress numbers (false)
%   thr		CEM2 threshold (1e-6)
%   animate	plot data and ellipses as the algorithm runs (false)
%   covtype	Covariance matrix structure: 1=diagonal, other=free (0)
% At least Cmax should be set explicitly.
% Example:
%   estS = gmmb_fj(data, 'Cmax', 50, 'thr', 1e-9)
%
% References:
%   [1] Duda, R.O., Hart, P.E, Stork, D.G, Pattern Classification,
%   2nd ed., John Wiley & Sons, Inc., 2001.
%   [2] Bilmes, J.A., A Gentle Tutorial of the EM Algorithm and its
%    Application to Parameter Estimation for Gaussian Mixture and Hidden
%    Markov Models
%   International Computer Science Institute, 1998
%   [3] Figueiredo, M.A.T., Jain, A.K., Unsupervised Learning on
%    Finite Mixture Models, IEEE transactions of pattern analysis and
%    machine intelligence, vol.24, no3, March 2002
%
% This code is directly based on [3] and code published on
% Figueiredo homepage: http://www.lx.it.pt/~mtf/
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
%   V_0_1 1.7  2003/09/17 18:20:19
%

function [estimate, varargout] = gmmb_fj(data, varargin);

[N, D] = size(data);	% number of points (n), dimensions (d)

% defaults
C = ceil(min(50, N/(D*D)/3));	% max n of components (k_max)
Cmin = 1;			% k_min
verbose = 0;
threshold = 1e-6;
animation = 0;
covtype = 0; % 1: diagonal, other: free

% parameter parsing
argc = nargin-1;
for argl = 1:2:(argc-1)
	switch lower(varargin{argl})
	  case 'cmax'
	    C = varargin{argl+1};
	  case 'cmin'
	    Cmin = varargin{argl+1};
	  case 'verbose'
	    verbose = varargin{argl+1};
	  case 'thr'
	    threshold = varargin{argl+1};
	  case 'animate'
	    animation = varargin{argl+1};
	  case 'covtype'
	    covtype = varargin{argl+1};
	  otherwise
	end
end

if C<1 ||C>N
	C = N;
	mu = data.';
else
	% initialize mu as random points from data
	permi = randperm(N);
	mu = data(permi(1:C),:).';  % D x C
end


% initialize sigma
s2 = max(diag(covfixer2(cov(data,1))/10));
sigma = repmat(s2*eye(D), [1 1 C]);

% weights initialization
alpha = ones(1,C) * (1/C);

if covtype == 1
	Nparc = D+D;	% (N)
else
	Nparc = D+D*(D+1)/2;	% (N)
end
Nparc2 = Nparc/2;

cost = repmat([-1], 1, 1000);

if animation
	aniH = my_plot_init;
	my_plot_ellipses(aniH, data, mu, sigma, alpha);
end


t = 0;
Cnz = C;	% (k_nz) k = kmax
Lmin = NaN;

u = zeros(N,C);	% semi_indic.'
for c = 1:C
	u(:,c) = cmvnpdf(data, mu(:,c).', sigma(:,:,c));
end
indic = u .* repmat(alpha, N,1);

old_loglike = sum(log(sum(realmin+indic, 2)));
old_L = Nparc2*sum(log(alpha)) + (Nparc2+0.5)*Cnz*log(N) - old_loglike;


while Cnz >= Cmin
	repeating = 1;
	
	while repeating
		t = t+1;
		
		c = 1;
		while c <= C
			indic = u .* repmat(alpha, N,1);
			normindic = indic ./ (realmin + repmat(sum(indic,2), 1,C));
			
			normf = 1/sum(normindic(:,c));
			aux = repmat(normindic(:,c), 1,D) .* data;
			
			nmu = normf * sum(aux,1);
			mu(:,c) = nmu.';
			
			if covtype == 1
				nsigma =  normf*diag(sum(aux .* conj(data), 1)) - diag(nmu.*conj(nmu));
			else
				nsigma =  normf*(aux' * data) - nmu'*nmu;
			end
			sigma(:,:,c) = covfixer2(nsigma);
			
			alpha(c) = max(0, sum(normindic(:,c))-Nparc2) / N;
			alpha = alpha / sum(alpha);
			

			if alpha(c) == 0
				Cnz = Cnz -1;
			else
				try
					u(:,c) = cmvnpdf(data, mu(:,c).', sigma(:,:,c));
				catch
					disp('covariance went bzrk !!!');
					sigma(:,:,c)
					%keyboard
					Cnz = 0;
				end
			end
			c=c+1;
			
			if Cnz == 0
				% number of components fell to zero
				% nothing can be done, return error
				estimate = NaN;
				return
			end

		end % while c <= C

		% purge alpha == 0 if necessary
		if length(find(alpha==0)) > 0
			nz = find(alpha>0);
			alpha = alpha(nz);
			mu = mu(:,nz);
			sigma = sigma(:,:,nz);
			u = u(:,nz);
			C = length(nz);
		end

		if animation
			my_plot_ellipses(aniH, data, mu, sigma, alpha);
		end
		
		u = zeros(N,C);	% semi_indic.'
		for c = 1:C
			u(:,c) = cmvnpdf(data, mu(:,c).', sigma(:,:,c));
		end
		indic = u .* repmat(alpha, N,1);
		
		loglike = sum(log(realmin+sum(indic, 2)));
		L = Nparc2*sum(log(alpha)) + (Nparc2+0.5)*Cnz*log(N) - loglike;

		
		if t<1001
			cost(t) = L;
		end

		if verbose
			disp(['Cnz=' num2str(Cnz) ' t=' num2str(t) ' ' num2str(abs(loglike - old_loglike)) ...
				' <? ' num2str(threshold*abs(old_loglike))]);
			disp(['t=' num2str(t) ' L= ' num2str(L)]);
		end

		diff = abs(loglike - old_loglike);
		if (diff < threshold*abs(old_loglike))
			repeating = 0;
		end
		
		old_L = L;
		old_loglike = loglike;
	end % while repeating
	
	if isnan(Lmin) || (L <= Lmin)
		Lmin = L;
		estimate = struct('mu', mu,...
			'sigma', sigma,...
			'weight', alpha.');
	end
	if verbose
		disp(['Cnz = ' num2str(Cnz)]);
	end

	% annihilate the least probable component
	m = find(alpha == min(alpha(alpha>0)));
	alpha(m(1)) = 0;
	Cnz = Cnz -1;
	% alpha doesn't need to be normalized here, even if it would seem logical to do so.
	
	if Cnz > 0
		alpha = alpha / sum(alpha);
	
		% purge alpha == 0 if necessary
		if length(find(alpha==0)) > 0
			nz = find(alpha>0);
			alpha = alpha(nz);
			mu = mu(:,nz);
			sigma = sigma(:,:,nz);
			u = u(:,nz);
			C = length(nz);
		end
		
		u = zeros(N,C);	% semi_indic.'
		for c = 1:C
			u(:,c) = cmvnpdf(data, mu(:,c).', sigma(:,:,c));
		end
		indic = u .* repmat(alpha, N,1);
		
		old_loglike = sum(log(realmin+sum(indic, 2)));
		old_L = Nparc2*sum(log(alpha)) + (Nparc2+0.5)*Cnz*log(N) - old_loglike;
	end
end

if(nargout>1)
	varargout(1) = {[Cnz t]};
end
if(nargout>2)
	varargout(2) = {cost(cost~=-1)};
end

% purge alpha==0
e = estimate;
inds = find(e.weight>0);
estimate.mu = e.mu(:,inds);
estimate.sigma = e.sigma(:,:,inds);
estimate.weight = e.weight(inds);

if animation
	my_plot_ellipses(aniH, data, estimate.mu, estimate.sigma, estimate.weight);
end

%disp(['Cfinal = ' num2str(length(inds))]);

% -----------------------------------------------------------

function h = my_plot_init;
h = figure;
figure(h);
title('Distribution of x_1 and x_2 values','FontSize',14);
xlabel('x_1 value','FontSize',14);
ylabel('x_2 value','FontSize',14);
zlabel('weight','FontSize',14);
view(2)
tic;

function my_plot_ellipses(h, data, mu, sigma, weight);
dtime = 0.3;

D = size(mu, 1);

if D ~= 2
	error('Can plot only 2D objects.');
end

[x,y,z] = cylinder([2 2], 40);
xy = [ x(1,:) ; y(1,:) ];

figure(h);

plot(data(:,1), data(:,2), 'rx');

hold on
C = size(mu, 2);
for c = 1:C
	mxy = chol(sigma(:,:,c))' * xy;
	x = mxy(1,:) + mu(1,c);
	y = mxy(2,:) + mu(2,c);
	z = ones(size(x))*weight(c);
	plot3(x,y,z, 'k-');
end
drawnow;
hold off

t = toc;
if t+0.01<dtime
	pause(dtime-t);
end
tic

