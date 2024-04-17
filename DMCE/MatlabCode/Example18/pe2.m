x=[0:0.1:1]';
x0=[0:0.0005:1]';
y0=2*x0.*cos(4*pi*x0);
y=2*x.*cos(4*pi*x);



n=length(x);
xx = ones(n,1);
b = mean(y);
 
pe(1) = mean((y0-mean(y)).^2);
xx0 = ones(size(x0));
for k=2:11
   xx=[xx,(x-0.5).^(k-1)];
   b =pinv(xx'*xx)*xx'*y;
   xx0 = [xx0,(x0-0.5).^(k-1)];
   fh = xx0*b ; 
 
	pe(k) = mean((fh-y0).^2);
   if k==10;
      bb0=b;
      xx00=xx;
      xx000=xx0;
   end;
   
end;

a=[0:10];
figure(2)
 
plot(a,pe,'-')


%xx = [ones(n,1),x,x.^2,x.^3,x.^5];
b = inv(xx00'*xx00)*xx00'*y;
yh = xx00*b;

%xx0 = [ones(size(x0)),x0,x0.^2,x0.^3,x0.^5];
yh0 = xx000*b; 

figure(1)

plot(x0,y0,'-',x,y,'o',x0,yh0,':',x,yh,'*')
legend('True', 'Sample','Fitted Curve','Fitted Values')