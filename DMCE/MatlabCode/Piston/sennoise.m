%p12.m

load('piston12.m')
data = piston12;
y=data(:,8);
x=data(:,2:7);
%[a,b,c,d]=krig(x,y)

[n,d]=size(x);

y = y + 0.5*std(y)*randn(n,1);
stdx=std(x)';
 

a=3.7;

lambda=0.5*sqrt(log(n)/n);
 
for k=1:60
   gamma = 0.01*1.1^k./stdx;
   ell0(k)=kriglkhd(x,y,gamma);
   p1 = lambda^2 - (abs(gamma) - lambda).^2.*(abs(gamma) < lambda);
   ell1(k) = ell0(k) - n*sum(p1);
end;
b0=0.01*1.1.^[1:60];

 
subplot(221)
plot(b0,ell1)
axis0=axis;
axis([0,3,axis0(3),axis0(4)]);
xlabel('\theta')
ylabel('Penalized Log-lkhd')
title('(a) N(0,0.25\sigma^2)')


y = y + 1*std(y)*randn(n,1);
stdx=std(x)';
 

a=3.7;

lambda=0.5*sqrt(log(n)/n);
 
for k=1:60
   gamma = 0.01*1.1^k./stdx;
   ell0(k)=kriglkhd(x,y,gamma);
   p1 = lambda^2 - (abs(gamma) - lambda).^2.*(abs(gamma) < lambda);
   ell1(k) = ell0(k) - n*sum(p1);
end;
b0=0.01*1.1.^[1:60];

 
subplot(222)
plot(b0,ell1)
axis0=axis;
axis([0,3,axis0(3),axis0(4)]);
xlabel('\theta')
ylabel('Penalized Log-lkhd')
title('(b) N(0,\sigma^2)')


y = y + sqrt(2)*std(y)*randn(n,1);
stdx=std(x)';
 

a=3.7;

lambda=0.5*sqrt(log(n)/n);
 
for k=1:60
   gamma = 0.01*1.1^k./stdx;
   ell0(k)=kriglkhd(x,y,gamma);
   p1 = lambda^2 - (abs(gamma) - lambda).^2.*(abs(gamma) < lambda);
   ell1(k) = ell0(k) - n*sum(p1);
end;
b0=0.01*1.1.^[1:60];

 
subplot(223)
plot(b0,ell1)
axis0=axis;
axis([0,3,axis0(3),axis0(4)]);
xlabel('\theta')
ylabel('Penalized Log-lkhd')
title('(c) N(0,2\sigma^2)')


y = y + sqrt(3)*std(y)*randn(n,1);
stdx=std(x)';
 

a=3.7;

lambda=0.5*sqrt(log(n)/n);
 
for k=1:60
   gamma = 0.01*1.1^k./stdx;
   ell0(k)=kriglkhd(x,y,gamma);
   p1 = lambda^2 - (abs(gamma) - lambda).^2.*(abs(gamma) < lambda);
   ell1(k) = ell0(k) - n*sum(p1);
end;
b0=0.01*1.1.^[1:60];

 
subplot(224)
plot(b0,ell1)
axis0=axis;
axis([0,3,axis0(3),axis0(4)]);
xlabel('\theta')
ylabel('Penalized Log-lkhd')
title('(d) N(0,3\sigma^2)')


