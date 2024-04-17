function [mu_est,ga_est,beta_est,sigma2_est, ell,step] = penkrigscad(x,y,lambda,ga0);
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
%			ga0 -- an initial value of gamma.
%
% Outputs:
%			mu_est = the estimate of mu
%			ga_est = the estimate of gamma
%			beta_est = the estiamte of beta used for kriging: y = hat{mu} + r(x,X)* beta
%			sigma2_est = the estimate of sigma^2
%			ell = the value of penalized likelihood function corresponding the final estimate
%			step = the number of iterations.
%

q=2;
epsr=1e-5; 

[n,d]=size(x);

x0 = zeros(n,n,d);

for k=1:d
   t1 = x(:,k);
   x0(:,:,k) = t1*ones(1,n)-ones(n,1)*t1';
end;
x0 = -abs(x0).^q;

delta = 2;
eps=1*10^(-3);
 

ell0=kriglkhd(x,y,ga0)-n*lambda*sum( ga0.^2 );

R00 = zeros(n,n);
 
for k=1:d
    R00 = R00 + ga0(k)*x0(:,:,k);
end;

R0 = exp(R00);
   
mu = (ones(1,n)*pinv(R0)*ones(n,1))^(-1)*ones(1,n)*pinv(R0)*y;
   
e = y-mu;
sigma2 = e'*pinv(R0)*e/n;


%% set the estimate corresponding to the initial values.

ell = ell0;
mu_est = mu;
sigma2_est = sigma2;
ga_est =ga0;
step =  0;

theta0 = log(ga0);

while (delta > eps) & (step<100)
   step=step+1;
   R01 = zeros(n,n,d);
   
   for k=1:d  % To compute partial derivative of R
      R01(:,:,k) = R0.*x0(:,:,k) *ga0(k);
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
   p1 = 2*n*lambda*ga0.^2;
	p2 = 4*n*lambda*diag(ga0.^2);
   
   if step <=10
       theta0 = theta0 + 0.5*pinv(q2 + p2,epsr)*(q1 - p1);
    else
       theta0 = theta0 + 0.25*pinv(q2 + p2,epsr)*(q1 - p1);
   end;
   
      
   ga0 =exp(theta0);
   
   ell1=kriglkhd(x,y,ga0)-n*lambda*sum( ga0.^2 );

	delta = abs((ell1-ell0)/ell0);
   ell0=ell1;
   
   R00 = zeros(n,n);
 
   for k=1:d
      R00 = R00 + ga0(k)*x0(:,:,k);
   end;
   R0 = exp(R00);
   mu = (ones(1,n)*pinv(R0)*ones(n,1))^(-1)*ones(1,n)*pinv(R0)*y;
   e = y-mu;
   sigma2 = e'*pinv(R0)*e/n;
   
   if ell0 > ell;
      ell = ell0;
      mu_est = mu;
      ga_est = ga0;
      sigma2_est = sigma2;
   end;
   
end;

 
beta_est = pinv(R0)*(y-mu_est);  
  
   

