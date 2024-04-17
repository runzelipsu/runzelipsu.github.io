function [lambda_sel,gcv]=gcv_scad(data,s0);
% [lambda_sel,gcv]=gcv_n(data) is used to selected
% the thresholding parameter lambda  for new thresholding with
% a=2+sqrt(3) considered in Fan and Li (1999). 
% Variable selection via Penalized likelihood.
%
% Input: data -- the data matrix arranged in the format [x,y];
%
% Output: lambda_sel -- the selected parameter lambda,
%          gcv -- the generalized cross-validation statistics
a=3.7;
[n,d]=size(data);
x=data(:,1:d-1);
y=data(:,d);
epsilon=10^(-4);
beta_int=inv(x'*x)*x'*y;
%l0=[0.25:0.5:10]*sqrt(s0)/sqrt(n);
l0=[0.05:0.15:10]*sqrt(s0)/sqrt(n);

gcv=zeros(1,length(l0));
for k=1:length(l0);
   lambda=l0(k);
    beta=scad(beta_int,lambda,data,s0);  
    rss=sum((y-x*beta).^2);
    if sum(beta~=0)>0,
       x0=x(:,(beta~=0));
       abeta=abs(beta(beta~=0));
       pbeta=(abeta<=lambda)+((a*lambda-abeta)/((a-1)*lambda)).*((abeta>lambda).*(abeta<a*lambda));
       pt=trace(x0*(x0'*x0+n*lambda*diag(pbeta./abeta))^(-1)*x0');
       gcv(1,k)=rss/(1-pt/n)/(1-pt/n)/n;
    else
       pt=d-1;
    end;   
   
end;
[med,medi]=min(gcv);
lambda_sel=l0(medi(1)); 

