function [mu_est,ga_est,step] = krig(x,y);
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
 
 
[n,d]=size(x);

x0 = zeros(n,n,d);

for k=1:d
   t1 = x(:,k);
   x0(:,:,k) = t1*ones(1,n)-ones(n,1)*t1';
end;
x0 = -abs(x0).^q;

mu0 = mean(y);
%ga0 = 0.5*1./range(x)';
ga0 =0.05*ones(d,1);
 

ga0(d+1)=var(y);

R00 = zeros(n,n);
 
 for k=1:d
    R00 = R00 + ga0(k)*x0(:,:,k);
 end;
 
R0 = ga0(d+1)*exp(R00);
 
mu = (ones(1,n)*pinv(R0)*ones(n,1))^(-1)*ones(1,n)*pinv(R0)*y;
e = y-mu;
  
ell0 = -0.5*log(det(R0)) -0.5*e'*pinv(R0)*e
  


delta = 2;
eps=10^(-4);
step = 0;

while (delta > eps) & (step<100)
   step=step+1;
    
   
   
   R01 = zeros(n,n,d+1);
   
   for k=1:d  % To compute partial derivative of R
      R01(:,:,k) = R00.*x0(:,:,k);
   end;
   
   R01(:,:,d+1) = R00;
   
   
   q1 = zeros(d+1,1);
   q2 = zeros(d+1,d+1);
   
   for k=1:d+1  % To compute the first partial derivatives
      q1(k) = 0.5*trace(pinv(R0)*(e*e' - R0)*pinv(R0)*R01(:,:,k));
   end;
   
   for k=1:d+1 % To compute the Hessian matrix
      for j=k:d+1
         q2(k,j) = 0.5*trace(pinv(R0)*R01(:,:,k)*pinv(R0)*R01(:,:,j));
         q2(j,k) = q2(k,j);
      end;
   end;
   
   ga0 = ga0 + pinv(q2)*q1;
   
   R00 = zeros(n,n);
 
   for k=1:d
      R00 = R00 + ga0(k)*x0(:,:,k);
   end;
   R0 = ga0(d+1)*exp(R00);
   
   mu = (ones(1,n)*pinv(R0)*ones(n,1))^(-1)*ones(1,n)*pinv(R0)*y;
   e = y-mu;
   %   sigma2 = e'*pinv(R0)*e/n;
	ell = -0.5*log(det(R0)) -0.5*e'*pinv(R0)*e

   
 
	delta =  max(abs(q1));
 
   
  end;

mu_est= mu;
%sigma2_est=sigma2;
ga_est = ga0;

