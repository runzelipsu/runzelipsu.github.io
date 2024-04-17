function [beta,std_beta,step]=lasso(beta_int,lambda,data,s0)
%inputs:
%       data=[x,y];
%       beta_int: initial value of beta, a d-1 column vector;
%
[n,d]=size(data);
x=data(:,1:d-1);
y=data(:,d);
epsilon=10^(-4);
delta=2;
step=0;
beta=zeros(d-1,1);
  
th0=2.5*sqrt( s0*diag(inv(x'*x)) ) ;
         % thresholding parameter for beta=0

while (delta>=epsilon)&(step<20)
 step=step+1;
 beta(abs(beta_int)<=th0)=zeros(sum(abs(beta_int)<=th0),1);
 x0=x(:,abs(beta_int)>th0);
 beta0=beta_int(abs(beta_int)>th0);  
 abeta=abs(beta0);
 pbeta0=ones(length(abeta),1);
 if length(beta0)>0
  beta0=(x0'*x0+n*lambda*diag(pbeta0./(abeta)))^(-1)*x0'*y;
  beta(abs(beta_int)>th0)=beta0;
  delta=max(abs(beta-beta_int));
  beta_int=beta;
 else
  delta=epsilon/2;
 end;

 % updating th0: keep the same for beta=0 and update the corresponding
 % th0 of beta~=0;

if sum(beta~=0)>0
 x0=x(:,(beta~=0));
 s0=sum((y-x0*inv(x0'*x0)*x0'*y).^2)/(n-sum(beta~=0)); % sigma^2
 th0(beta~=0)=1.5*sqrt( s0*diag(inv(x0'*x0)) ) ;
end;
end;


std_beta=zeros(d-1,1);

if sum(beta~=0)>0
   s0=sum((x*beta-y).^2)/(n-d+1);
   x0=x(:,(beta~=0));
   abeta=abs(beta0);
   pbeta0=ones(length(abeta),1);
 
   temp=diag((x0'*x0+n*lambda*diag(pbeta0./abeta))^(-1)...
                           *x0'*x0*(x0'*x0+n*lambda*diag(pbeta0./abeta))^(-1));
   std_beta(beta~=0)=sqrt(s0*temp);
end;





 
