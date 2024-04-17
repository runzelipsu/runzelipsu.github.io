function [mu_est,ga_est,beta_est,Psi_est, ell,step]  = mvkrigmle(x,y,ga0);
%
% [mu_est, ga_est] = krigmle(x,y,q) is to compute the MLE of
% mu and gamma in Li and Sudjianto (2003) Penalized likelihood for 
% Gaussian kriging models with applications to analysis of computer
% experiments
%
% Inputs: 
%			x -- x variables
%			y --  response variable
%			ga0 -- an initial value of gamma.
%
% Outputs:
%			mu_est = the estimate of mu
%			ga_est = the estimate of gamma
%

q=2;
epsr=1e-5; 

[n,d]=size(x);
[r,s]=size(y);

x0 = zeros(n,n,d);

for k=1:d
   t1 = x(:,k);
   x0(:,:,k) = t1*ones(1,n)-ones(n,1)*t1';
end;
x0 = -abs(x0).^q;

 
delta = 2;
eps=1*10^(-5);
 
ell0=mvkriglkhd(x,y,ga0);
 

R00 = zeros(n,n);
 
for k=1:d
    R00 = R00 + ga0(k)*x0(:,:,k);
end;

R0 = exp(R00);
   
mu = (ones(1,n)*pinv(R0)*ones(n,1))^(-1)*ones(1,n)*pinv(R0)*y;
   
e = y-ones(n,1)*mu;
Psi = e'*pinv(R0)*e/n;

%% set the estimate corresponding to the initial values.

ell = ell0;
mu_est = mu;
Psi_est = Psi;
ga_est =ga0;
step =  0;

theta0 = log(ga0);

theta_ini=theta0;

while (delta > eps) & (step<100)
   step=step+1;
    
   R01 = zeros(n,n,d);
   
   for k=1:d  % To compute partial derivative of R
      R01(:,:,k) = R0.*x0(:,:,k)*ga0(k);
   end;
   
 
   
   q1 = zeros(d,1);
   q2 = zeros(d,d);
   
   for k=1:d  % To compute the first partial derivatives
      q1(k) = trace(pinv(R0)*(e*pinv(Psi)*e' - s*R0)*pinv(R0)*R01(:,:,k));
   end;
   
   for k=1:d % To compute the Hessian matrix 
      for j=k:d
         q2(k,j) = -s*trace(pinv(R0,epsr)*R01(:,:,k)*pinv(R0,epsr)*R01(:,:,j));
         q2(j,k) = q2(k,j);
      end;
   end;
   
%   if step <=20
      theta0 = theta0 - pinv(q2,epsr)*q1;
%   else % step-halving
%      theta0 = theta0 - 0.5*pinv(q2,epsr)*q1;
%   end;
   
   
 	ga0 =exp(theta0);
   
   ell1 = mvkriglkhd(x,y,ga0);
   
   delta = abs((ell1-ell0)/ell0)
   
   ell0= ell1;   
   
    
   R00 = zeros(n,n);
 
   for k=1:d
      R00 = R00 + ga0(k)*x0(:,:,k);
   end;
   R0 = exp(R00);
   
   mu = (ones(1,n)*pinv(R0)*ones(n,1))^(-1)*ones(1,n)*pinv(R0)*y;
   
   e = y-ones(n,1)*mu;
   Psi = e'*pinv(R0)*e/n;
   
   if ell1 > ell;
      ell = ell1;
      mu_est = mu;
      ga_est = ga0;
      Psi_est = Psi;
   end;
   
  end;
  beta_est = pinv(R0)*(y-ones(r,1)*mu_est);  