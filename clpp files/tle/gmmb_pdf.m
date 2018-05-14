%GMMB_PDF    - (Complex range) multivariate Gaussian mixture model pdf
%
% p = gmmb_pdf(x, mu, sigma, weight)
%
% x = N x D vector
% mu = D x C matrix
% sigma = D x D x C matrix array
% weight = C x 1 vector
% p = N x 1 vector
%
% D dimensions, N points, C components

% Author: Pekka Paalanen <pekka.paalanen@lut.fi>

%
% V_0_1
% gmmb_pdf.m,v 1.4 2004/02/25 13:21:56 paalanen Exp

function [pdf] = gmmb_pdf(x, mu, sigma, weight);
N = size(x,1);
pdf = zeros(N,1);
for l = 1:size(weight,1);
	pdf = pdf + weight(l) * cmvnpdf(x, mu(:,l).', sigma(:,:,l));
end
