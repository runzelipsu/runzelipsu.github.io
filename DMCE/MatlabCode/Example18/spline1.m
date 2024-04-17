x=[0:0.1:1]';
x0=[0:0.01:1]';
y0=2*x0.*cos(4*pi*x0);
y=2*x.*cos(4*pi*x);
%figure(1)
%plot(x0,y0,'-',x,y,'o')

n=length(x);
xx = ones(n,1);
 

b1=[];
for i=1:7
   b1=[b1,(x-i/8).^2.*(x>i/8)];
end; 
 

xx = [ones(n,1),x,x.^2,b1];
b = inv(xx'*xx)*xx'*y;
yh = xx*b;

b2=[];

for i=1:7
   b2=[b2,(x0-i/8).^2.*(x0>i/8)];
end;

xx0 = [ones(size(x0)),x0,x0.^2,b2];
yh0 = xx0*b; 

figure(1)

plot(x0,y0,'-',x,y,'o',x0,yh0,':',x,yh,'*')
legend('True', 'Sample','Fitted Curve','Fitted Values')
xlabel('x')
ylabel('y')

figure(2)

plot(x0,b2,'-')