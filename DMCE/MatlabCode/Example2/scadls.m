function [beta,std_est,gcv,gcva,gcvb,effno,effnoa,effnob]=scadls(x,y,lambda, bi);

% x: the standardized design matrix
% y: the response vector
% bi: the initial value, LSE

a=3.7;

delta=1;

[n,d]=size(x);

v = x'*x/(n-1);
z=bi;

eps=10^(-4);
beta0=zeros(d,1);

step=0;
beta0=bi;

eff=zeros(d,0);
while (delta>eps) & (step<=100)
   
   step=step+1;
   
   for j=1:d
      zj =  beta0(j) - v(j,:) * (beta0 - z) ;
 
      if abs(zj) < lambda
         beta0(j) = 0;
         eff0(j)=0;
      end;
      if (lambda <= abs(zj) ) & (abs(zj) < a*lambda)
         beta0(j) = zj - sign(zj)*(a*lambda-abs(zj))/(a-1);
         eff0(j) = abs(beta0(j))/abs(z(j));
      end;
      if abs(zj)>= a*lambda
         beta0(j) = zj;
         eff0(j)=1;
      end;               
   end;
   effnoa = sum(eff0);
   delta = max(abs(beta0-bi));
   bi=beta0;
end;
beta=beta0;
effnob = sum(abs(beta)./abs(z));

betan = beta(beta~=0);
xn = x(:,beta~=0);
abeta=abs(betan);
p2 =-1/(a-1)*((abeta>=lambda).*(abeta<a*lambda));
Sigma=diag(p2);

q2I = inv(xn'*xn + n*Sigma);

effno = trace(q2I*(xn'*xn));
rss = mean((y-x*beta).^2);
gcv = rss/(1-effno/n)/(1-effno/n);
gcva = rss/(1-effnoa/n)/(1-effnoa/n);
gcvb = rss/(1-effnob/n)/(1-effnob/n);



sigma2 = rss*n/(n-effno);
q1 = xn'*xn;

std_est0 = sqrt(sigma2*diag(q2I*q1*q2I));
std_est=zeros(d,1);
std_est(beta~=0)=std_est0;

