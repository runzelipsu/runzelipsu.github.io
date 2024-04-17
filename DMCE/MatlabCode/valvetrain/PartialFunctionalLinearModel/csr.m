function [yhat, beta]=csr(x,y,n,xknot,p)

x2=x.*x;
x3=x2.*x;
xa=[ones(n,1), x, x2, x3]

for i=1:p
xknotv=(x-xknot(i)).*(x-xknot(i)).*(x-xknot(i));
xknotv=(xknotv+abs(xknotv))*0.5;
xa=[xa, xknotv];
end;

beta=inv(xa' * xa)*xa'*y;
yhat=xa*beta;
