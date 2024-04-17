function [phat, pci] = betafit(x,alpha)
%BETAFIT Parameter estimates and confidence intervals for beta distributed data.
%   BETAFIT(X) Returns the maximum likelihood estimates of the  
%   parameters of the beta distribution given the data in the vector, X.  
%
%   [PHAT, PCI] = BETAFIT(X,ALPHA) gives MLEs and 100(1-ALPHA) 
%   percent confidence intervals given the data. By default, the
%   optional paramter ALPHA = 0.05 corresponding to 95% confidence intervals.

%   Reference:
%      [1]  Hahn, Gerald J., & Shapiro, Samuel, S.
%      "Statistical Models in Engineering", Wiley Classics Library
%      John Wiley & Sons, New York. 1994 p. 95.

%   B.A. Jones 1-27-95
%   Copyright (c) 1993-98 by The MathWorks, Inc.
%   $Revision: 2.5 $  $Date: 1997/11/29 01:44:50 $

if nargin < 2 
    alpha = 0.05;
end
p_int = [alpha/2; 1-alpha/2];

if min(size(x)) > 1
  error('Requires a vector first input.');
end
 
n = length(x);

% Initial Estimates.
tmp1 = prod((1-x) .^ (1./n));
tmp2 = prod(x .^ (1./n));
tmp3 = (1 - tmp1 - tmp2);
ahat = 0.5*(1-tmp1) ./ tmp3;
bhat = 0.5*(1-tmp2) ./ tmp3;

pstart = [ahat bhat];
phat = fmins('betalike',pstart,[],[],x);

[logL,info]=betalike(phat,x);
sigma = sqrt(diag(info));
pci = norminv([p_int p_int],[phat; phat],[sigma';sigma']);
