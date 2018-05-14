function [evClass, evVal, varargout] = gmmb(classifierS, data, use_mahalanobis);
%GMMB - Run GMMBayes project classifier
%
% [evClass, evVal, varargout] = gmmb(classifierS, data, use_mahalanobis)
%   classifierS      classifier struct created by newgmmb()
%   data             samples-by-dimensions matrix of data points
%   use_mahalanobis  boolean flag, use mahalanobis distance or Bayesian
%                    classifier
%
%   evClass        vector of class numbers for classified data
%   evVal          the probability by which the class was selected (Bayesian)
%                  or negative mahalanobis distance from class center
%   varargout{1}   the full a posteriori probability matrix or
%                  mahalanobis distance matrix, samples-by-classes
%
% See also NEWGMMB, GMMB_CLASSIFY
%
% Author: Pekka Paalanen <pekka.paalanen@lut.fi>
%
%
% V_0_1
% gmmb.m,v 1.2 2004/02/25 13:21:56 paalanen Exp

if use_mahalanobis == 0
	P = gmmb_classify(classifierS.bayesS, data, 'values');
	% bigger probability is better
	[evVal, evClass] = max(P, [], 2);
else
	P = gmmb_classify(classifierS.bayesS, data, 'values', 'mahalanobis');
	% smaller distance is better
	[evVal, evClass] = min(P, [], 2);
	evVal = -evVal; % ...so reverse the order
end

if nargout>2
	varargout{1} = P;
end

