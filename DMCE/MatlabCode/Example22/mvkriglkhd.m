function [ell,mu,Psi] = mvkriglkhd(x,y,ga0);
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
[r,s]=size(y);

x0 = zeros(n,n,d);

for k=1:d
   t1 = x(:,k);
   x0(:,:,k) = t1*ones(1,n)-ones(n,1)*t1';
end;
x0 = -abs(x0).^q;

 
R00 = zeros(n,n);
 
 for k=1:d
    R00 = R00 + x0(:,:,k)*ga0(k);
 end;
 
 
R0 = exp(R00);

mu = (ones(1,n)*pinv(R0)*ones(n,1))^(-1)*ones(1,n)*pinv(R0)*y;
e = y-ones(n,1)*mu; % nxJ matrix

Psi = e'*pinv(R0)*e/n;

%R0=Sigma2*R0;

ell = -0.5*s*log(det(R0)) -0.5*r*log(det(Psi)) -0.5*trace(pinv(R0)*e*pinv(Psi)*e');

 


 