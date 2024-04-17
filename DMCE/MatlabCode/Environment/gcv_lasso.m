function [lambda_sel,gcv]=gcv_lasso(data,s0);
% [lambda_sel,gcv]=gcv_l(data) is used to selected
% the thresholding parameter lambda for lasso
% considered in Fan and Li (1999). Variable selection via Penalized
% likelihood.
%
% Input: data -- the data matrix arranged in the format [x,y];
%
% Output: lambda_sel -- the selected parameter lambda,
%                       for lasso,
%
%                gcv -- the generalized cross-validation statistics
[n,d]=size(data);
x=data(:,1:d-1);
y=data(:,d);
epsilon=10^(-4);
beta_int=inv(x'*x)*x'*y;
l0=[0.0005:0.1:2]*sqrt(s0)/sqrt(n);

for k=1:length(l0);
   lambda=l0(k);
   beta=lasso(beta_int,lambda,data,s0);
   %based on soft-threshold, coinciding with lasso;
   rss=sum((y-x*beta).^2);
   if sum(beta~=0)>0,
      x0=x(:,(beta~=0));
      abeta=abs(beta(beta~=0));
      pbeta=ones(length(abeta),1);
      pt=trace(x0*(x0'*x0+n*lambda*diag(pbeta./abeta))^(-1)*x0');
   else
      pt=d-1;
   end;   
   gcv(k,1)=rss/(1-pt/n)/(1-pt/n)/n;
 end;
[med,medi]=min(gcv);
lambda_sel=l0(medi(1));
 

