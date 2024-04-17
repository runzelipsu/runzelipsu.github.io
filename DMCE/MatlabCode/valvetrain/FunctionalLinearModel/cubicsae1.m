n=50;
p=10;
nl=15; 
dx=0.005;
displace=0.001;
randn('seed',1001);
rand('seed',1001)

xmu=0.0;
sigma2=1;
sigma=sqrt(sigma2);

dknot=1.0/(p+1.0);
xknot=[dknot:dknot:p*dknot]';

eps=sigma*randn(n,1)+xmu;
x=sort(rand(n,1));
y=4*sin(2*pi*x)+eps;;
%y=func2(x,eps);

%lambda0=logspace(log(0.001),log(1.5),nl)*sigma*log(n);
lambda0=logspace(log(0.001),log(1),nl)*sigma*log(n);
 
 
gcvv=[];
for i=1:nl
xx=x;
yy=y;
xlambda=lambda0(i);
[yhat, beta,gcv]=csr1(xx,yy,n,xknot,p,xlambda);
gcvv(i)=gcv;
inde(i)=xlambda;
end;

[gcvmin,indmin]=min(gcvv);
[yhat, beta,gcv]=csr1(xx,yy,n,xknot,p,inde(indmin));

[yhat1, beta1]=csr1(xx,yy,n,xknot,p,0);



xp=[0:0.001:1]';
yp=4*sin(2*pi*xp);
%yp=4*(xp-0.5).*(xp-0.5);

figure(1)
subplot(221)
plot(x, yhat,'-',xp,yp,':',x,yhat1,'--');
xlabel('t')
ylabel('Estimated \mu(t)')
title('(a) K=10')
%axis([0 1 -1 1]);
%grid on;
%print -dpsc csm1s4.eps;
subplot(222)
plot(log(lambda0),gcvv);
xlabel('log(\lambda)')
ylabel('GCV scores')
title('(b) K=10')
axis0=axis;
axis([axis0(1),0,axis0(3),axis0(4)])