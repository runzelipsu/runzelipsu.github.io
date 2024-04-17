%test2.m

load('piston12.m')
data = piston12;
y=data(:,8);
x=data(:,2:7);
%[a,b,c,d]=krig(x,y)

[n,d]=size(x);
stdx=std(x)';

a=3.7;

lambda=0.5*sqrt(log(n)/n);
lambda2=0.5*sqrt(log(n)/n);
lambda3=0.5*sqrt(log(n)/n);

for k=1:60
   gamma = 0.01*1.1^k./stdx;
   ell0(k)=kriglkhd(x,y,gamma);
   p1 = lambda^2 - (abs(gamma) - lambda).^2.*(abs(gamma) < lambda);
   ell1(k) = ell0(k) - n*sum(p1);
   p2=lambda2*abs(gamma);
   ell2(k) = ell0(k) - n *sum(p2);
   p3 = lambda3 * gamma.^2;
   ell3(k) = ell0(k) - n *sum(p3);
end;
b0=0.01*1.1.^[1:60];

subplot(221)
plot(b0,ell0)
axis0=axis;
axis([0,3,axis0(3),axis0(4)]);

subplot(222)
[tt1,tt2]=max(ell1);
b0(tt2)
plot(b0,ell1)
axis0=axis;
axis([0,3,axis0(3),axis0(4)]);


subplot(223)
[tt1,tt2]=max(ell2);
b0(tt2)

plot(b0,ell2)
axis0=axis;
axis([0,3,axis0(3),axis0(4)]);


subplot(224)
[tt1,tt2]=max(ell3);
b0(tt2)

plot(b0,ell3)
axis0=axis;
axis([0,3,axis0(3),axis0(4)]);

pause

[a1,b1,c1,d1]=penkrig(x,y,lambda)

gamma= c1;

ell50=kriglkhd(x,y,gamma);

p1 = lambda^2 - (abs(gamma) - lambda).^2.*(abs(gamma) < lambda);
ell51 = ell50 - n*sum(p1);

%[a1,b1,c1,d1]=penkrig(x,y,2)
