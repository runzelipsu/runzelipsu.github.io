function [mu_est,sigma2_est,ga_est,step] = krig(x,y);
%
% [mu_est, ga_est] = penkrig(x,y,q) is to compute the estimate of
% mu and gamma in Li and Sudjianto (2003) Penalized likelihood for 
% Gaussian kriging models with applications to analysis of computer
% experiments
%
% Inputs: 
%			x -- x variables
%			y --  response variable
%
% Outputs:
%			mu_est = the estimate of mu
%			ga_est = the estimate of gamma
%

q=2;
 
epsr =10^(-10);

[n,d]=size(x);

x0 = zeros(n,n,d);

for k=1:d
   t1 = x(:,k);
   x0(:,:,k) = t1*ones(1,n)-ones(n,1)*t1';
end;
x0 = -abs(x0).^q;

mu0 = mean(y);
ga0 = 1./(max(x)-min(x))';

delta = 2;
eps=10^(-2);
step = 0;

while (delta > eps) & (step<300)
   step=step+1
    
   R00 = zeros(n,n);
 
   for k=1:d
      R00 = R00 + ga0(k)*x0(:,:,k);
   end;
   R0 = exp(R00);
   
   mu = (ones(1,n)*pinv(R0)*ones(n,1))^(-1)*ones(1,n)*pinv(R0)*y;
   e = y-mu;
   sigma2 = e'*pinv(R0)*e/n;
   
   R01 = zeros(n,n,d);
   
   for k=1:d  % To compute partial derivative of R
      R01(:,:,k) = R00.*x0(:,:,k);
   end;
   
   q1 = zeros(d,1);
   q2 = zeros(d,d);
   
   for k=1:d  % To compute the first partial derivatives
      q1(k) = trace(pinv(R0)*(e*e'/sigma2 - R0)*pinv(R0)*R01(:,:,k));
   end;
   
   for k=1:d % To compute the Hessian matrix
      for j=k:d
         q2(k,j) = trace(pinv(R0)*R01(:,:,k)*pinv(R0)*R01(:,:,j));
         q2(j,k) = q2(k,j);
      end;
   end;
   
   
   if step <=20
      ga = ga0 + pinv(q2)*q1;
   else
      ga = ga0 + 0.5*pinv(q2)*q1;
   end;
   
      
   ga = ga.*(ga>=0);
	delta = max(abs(ga-ga0));
   ga0=ga;
end;

mu_est= mu;
sigma2_est=sigma2;
ga_est = ga;

