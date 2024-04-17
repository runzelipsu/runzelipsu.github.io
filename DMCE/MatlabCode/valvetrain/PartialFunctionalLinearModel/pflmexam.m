rand('seed',1001);
randn('seed',1001);
n=50;
J=80;
z=randn(n,J);
x=randn(n,3);
t=[1:J]/(J+1);
b0=ones(1,J);
b1=ones(1,J);
b2=sin(2*pi*t);
b3=4*(t-0.5).^2;
y = ones(n,1)*b0 +x(:,1)*b1+ x(:,2)*b2+x(:,3)*b3+z; % n x J matrix

xo=[ones(n,1),x];
beta_hat=inv(xo'*xo)*xo'*y;

p=10;
uknot=[1:p]'/(p+1);
ni=n;
nj=J;
nl=15;

sigma=1;

lambda=logspace(log(0.001),log(1.5),nl)*log(nj)*sigma/sqrt(ni);
%lambda2=logspace(log(0.001),log(1.5),nl)*log(nj)*sigma/sqrt(ni);

 
uj=t';

for i=1:nl
   [betas, beta,gcv]=psr(uj,beta_hat(1,:)',nj,uknot,p,lambda(i));
   gcvv0(i)=gcv;
   [betas, beta,gcv]=psr(uj,beta_hat(2,:)',nj,uknot,p,lambda(i));
   gcvv1(i)=gcv; 
   [betas, beta,gcv]=psr(uj,beta_hat(3,:)',nj,uknot,p,lambda(i));
   gcvv2(i)=gcv; 
   [betas, beta,gcv]=psr(uj,beta_hat(4,:)',nj,uknot,p,lambda(i));
   gcvv3(i)=gcv;    
end;

[gcvmin0,indmin0]=min(gcvv0);
[beta0s, beta,gcv]=psr(uj,beta_hat(1,:)',nj,uknot,p,lambda(indmin0));
[gcvmin1,indmin1]=min(gcvv1);
[beta1s, beta,gcv]=psr(uj,beta_hat(2,:)',nj,uknot,p,lambda(indmin1));
[gcvmin2,indmin2]=min(gcvv2);
[beta2s, beta,gcv]=psr(uj,beta_hat(3,:)',nj,uknot,p,lambda(indmin2));
[gcvmin1,indmin3]=min(gcvv3);
[beta3s, beta,gcv]=psr(uj,beta_hat(4,:)',nj,uknot,p,lambda(indmin3));


delta=2;
eps=5*10^(-5);
step=0;
x1=xo(:,1:2);
x2=xo(:,3:4);
b_ini=[2,2]';
while (delta>eps)&(step<=200)
   step=step+1;
   y1 = y - x2(:,1)*beta2s' - x2(:,2)*beta3s';
   y1m = (mean(y1'))';
   b=inv(x1'*x1)*x1'*y1m;
   y2=y-(x1*b)*ones(1,J);
   beta_h=inv(x2'*x2)*x2'*y2;
	for i=1:nl
   	[betas, beta,gcv]=psr(uj,beta_h(1,:)',nj,uknot,p,lambda(i));
	   gcvv0(i)=gcv;
   	[betas, beta,gcv]=psr(uj,beta_h(2,:)',nj,uknot,p,lambda(i));
	   gcvv1(i)=gcv; 
	end;

	[gcvmin0,indmin0]=min(gcvv0);
	[beta2s, beta,gcv]=psr(uj,beta_h(1,:)',nj,uknot,p,lambda(indmin0));
	[gcvmin1,indmin1]=min(gcvv1);
	[beta3s, beta,gcv]=psr(uj,beta_h(2,:)',nj,uknot,p,lambda(indmin1));
   delta=max(abs(b-b_ini));
   b_ini=b;
end;

b'

figure(1)
subplot(2,1,1)
plot(uj,beta2s,'-',t,b2,':')
xlabel('t')
ylabel('\beta_2(t)')
title('(a) Estimate of \beta_2(t)')

subplot(2,1,2)
plot(uj,beta3s,'-',t,b3,':')
xlabel('t')
ylabel('\beta_3(t)')
title('(b) Estimate of \beta_3(t)')


   
