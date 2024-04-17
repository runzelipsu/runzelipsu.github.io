function [mu_est,sigma2_est,ga_est,step] = penkrig(x,y,lambda);
%
% [mu_est, ga_est] = penkrig(x,y,q) is to compute the estimate of
% mu and gamma in Li and Sudjianto (2003) Penalized likelihood for 
% Gaussian kriging models with applications to analysis of computer
% experiments
%
% Inputs: 
%			x -- x variables
%			y --  response variable
%			lambda --  a tuning parameter;
%
% Outputs:
%			mu_est = the estimate of mu
%			ga_est = the estimate of gamma
%

q=2;
a=3.7;
epsr =0;

[n,d]=size(x);

x0 = zeros(n,n,d);

for k=1:d
   t1 = x(:,k);
   x0(:,:,k) = t1*ones(1,n)-ones(n,1)*t1';
end;
x0 = -abs(x0).^q;

 
ga0 = 0.01*1.1^15./std(x)';

delta = 2;
eps=10^(-4);
step = 0;

ell0=kriglkhd(x,y,ga0) ...
   -n*sum( lambda^2 - (abs(ga0) - lambda).^2.*(abs(ga0) < lambda));


while (delta > eps) & (step<100)
   step=step+1;
    
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
   
   aga =abs(ga0);
   p10 =(aga<=lambda)+((a*lambda-aga)/((a-1)*lambda)).*((aga>lambda).*(aga<a*lambda));
   p1 = n*lambda*((aga.*p10)./(aga+epsr));
	p2 = n*lambda*diag(aga./(aga+epsr));
   ga0 = ga0 + pinv(q2 + p2)*(q1 - p1);
   
   ell1=kriglkhd(x,y,ga0) ...
	   -n*sum( lambda^2 - (abs(ga0) - lambda).^2.*(abs(ga0) < lambda));


	delta = abs((ell1-ell0)/ell0)
   if ell1>ell0;
      ell0 = ell1;
   end;
   
end;

mu_est= mu;
sigma2_est=sigma2;
ga_est = ga0;

