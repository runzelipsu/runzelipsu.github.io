%p12.m

load('piston12.m')
data = piston12;
y=data(:,8);
x=data(:,2:7);
%[a,b,c,d]=krig(x,y)

[n,d]=size(x);

y = y + 0*std(y)*randn(n,1);
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
title('(a) \lambda = \lambda_1')


lambda=0.75*sqrt(log(n)/n);
 
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
title('(b) \lambda = \lambda_2')



lambda=1*sqrt(log(n)/n);
 
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
title('(c) \lambda = \lambda_3')



lambda=1.25*sqrt(log(n)/n);
 
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
title('(b) \lambda = \lambda_4')




 
