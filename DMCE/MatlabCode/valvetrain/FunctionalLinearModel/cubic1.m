n=200;
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
y=func1(x,eps);
%y=func2(x,eps);

%lambda0=logspace(log(0.001),log(1.5),nl)*sigma*log(n);
lambda0=logspace(log(0.001),log(1.5),nl)*sigma*log(n);
 
 
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


xp=[0:0.001:1]';
yp=4*sin(2.0*pi*xp);
%yp=4*(xp-0.5).*(xp-0.5);

figure(1)
plot(x, yhat,'-',xp,yp,':');
axis([0 1 -4 4]);
%grid on;
%print -dpsc csm1s4.eps;
figure(2)
plot(gcvv);
