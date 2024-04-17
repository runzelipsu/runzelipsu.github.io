function [ypen, beta,gcv]=psr(x,y,n,xknot,p,xlambda)

x2=x.*x;
x3=x2.*x;
xa=[ones(n,1), x, x2, x3];

sigmam=zeros(p+4,p+4);
sigmam(5:p+4,5:p+4) = xlambda*eye(p);

n=length(x);
p=length(xknot);

xknotv=(x*ones(1,p)-ones(n,1)*xknot').^3;
xknotv=(xknotv+abs(xknotv))*0.5;

xa=[xa, xknotv];
 
h=inv(xa' * xa+sigmam)*xa';
beta=h*y;
h=xa*h;
ypen=h*y;

elambda=trace(h);
gcv=n*sum((y-ypen).^2)/(n-elambda)/(n-elambda);
