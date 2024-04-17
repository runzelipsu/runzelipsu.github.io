rand('seed',0);
randn('seed',0);
n=200;
p=10;
nl=10;
displace=0.001;
xmu=0.0;

sigma=4;
sigma2=sigma*sigma;

dknot=1.0/(p+1.0);
xknot=[dknot:dknot:p*dknot]';

eps=sigma*randn(n,1)+xmu;

x=sort(rand(n,1));
%y=func1(x,eps);
y=func2(x,eps);

%Calculate the cubic spline regression and store it to ycsr
[ycsr, beta]=csr(x,y,n,xknot,p);

%Calculate the penalized regression spline and store it to ypen
lambda0=logspace(log(0.001),log(1.5),nl)*sigma*log(n);

gcvv=[];
for i=1:nl
xx=x;
yy=y;
xlambda=lambda0(i);
[ypen, beta,gcv]=psr(xx,yy,n,xknot,p,xlambda);
gcvv(i)=gcv;
inde(i)=xlambda;
end;

[gcvmin,indmin]=min(gcvv);
[ypen, beta,gcv]=psr(xx,yy,n,xknot,p,inde(indmin));

xp=[0:0.001:1]';
%yp=4*sin(2.0*pi*xp);
yp=16*(xp-0.5).*(xp-0.5);

figure(1);
plot(xp,yp,'k:',x, ycsr,'k-.',x,ypen,'k-');
legend('True Function','Cubic Spline','Penalized Spline',0);
axis([0 1 -0.5 4.5]);
%grid on;
print -dpsc csm2s4.eps;

figure(2);
plot(gcvv)


p=5;
dknot=1.0/(p+1.0);
xknot=[dknot:dknot:p*dknot]';
for i=1:nl
xx=x;
yy=y;
xlambda=lambda0(i);
[y5, beta,gcv]=psr(xx,yy,n,xknot,p,xlambda);
gcvv(i)=gcv;
inde(i)=xlambda;
end;

[gcvmin,indmin]=min(gcvv);
[y5, beta,gcv]=psr(xx,yy,n,xknot,p,inde(indmin));

p=10;
dknot=1.0/(p+1.0);
xknot=[dknot:dknot:p*dknot]';
for i=1:nl
xx=x;
yy=y;
xlambda=lambda0(i);
[y10, beta,gcv]=psr(xx,yy,n,xknot,p,xlambda);
gcvv(i)=gcv;
inde(i)=xlambda;
end;

[gcvmin,indmin]=min(gcvv);
[y10, beta,gcv]=psr(xx,yy,n,xknot,p,inde(indmin));

p=15;
dknot=1.0/(p+1.0);
xknot=[dknot:dknot:p*dknot]';
for i=1:nl
xx=x;
yy=y;
xlambda=lambda0(i);
[y15, beta,gcv]=psr(xx,yy,n,xknot,p,xlambda);
gcvv(i)=gcv;
inde(i)=xlambda;
end;

[gcvmin,indmin]=min(gcvv);
[y15, beta,gcv]=psr(xx,yy,n,xknot,p,inde(indmin));

p=20;
dknot=1.0/(p+1.0);
xknot=[dknot:dknot:p*dknot]';
for i=1:nl
xx=x;
yy=y;
xlambda=lambda0(i);
[y20, beta,gcv]=psr(xx,yy,n,xknot,p,xlambda);
gcvv(i)=gcv;
inde(i)=xlambda;
end;

[gcvmin,indmin]=min(gcvv);
[y20, beta,gcv]=psr(xx,yy,n,xknot,p,inde(indmin));

figure(3);
plot(xp,yp,'k:',x, y5,'k-',x,y10,'k-.',x,y15,'k--',x,y20,'k.');
%plot(xp,yp,'k-',x, y5,'k.',x,y10,'r.',x,y15,'c.');
legend('True Function','Penalized Spline K=5','Penalized Spline K=10','Penalized Spline K=15','Penalized Spline K=20',0);
%axis([0 1 -4.5 4.5]);
%print -dpsc csm1s4_kcompa.eps;

