%p12.m

load('piston12.m')
data = piston12;
y=data(:,8);
x=data(:,2:7);
 

[n,d]=size(x);

 
stdx=std(x)';
 

a=3.7;

lambda=0.5*sqrt(log(n)/n);
lambda2=0.5*sqrt(log(n)/n);
lambda3=0.5*sqrt(log(n)/n);

for k=1:120
   gamma = 0.01*1.05^k./stdx;
   ell0(k)=kriglkhd(x,y,gamma);
   p1 = lambda^2 - (abs(gamma) - lambda).^2.*(abs(gamma) < lambda);
   ell1(k) = ell0(k) - n*sum(p1);
   p2=lambda2*abs(gamma);
   ell2(k) = ell0(k) - n *sum(p2);
   p3 = lambda3 * gamma.^2;
   ell3(k) = ell0(k) - n *sum(p3);
end;
b0=0.01*1.05.^[1:120];

subplot(221)
plot(b0,ell0)
axis0=axis;
axis([0,3,axis0(3),axis0(4)]);
xlabel('\theta')
ylabel('Log-likelihood')
title('(a) Log-likelihood Function')

subplot(222)
plot(b0,ell1)
axis0=axis;
axis([0,3,axis0(3),axis0(4)]);
xlabel('\theta')
ylabel('Penalized Log-lkhd')
title('(b) SCAD Penalized Log-likelihood')



subplot(223)
plot(b0,ell2)
axis0=axis;
axis([0,3,axis0(3),axis0(4)]);
xlabel('\theta')
ylabel('Penalized Log-lkhd')
title('(c) L_1 Penalized Log-likelihood')



subplot(224)
plot(b0,ell3)
axis0=axis;
axis([0,3,axis0(3),axis0(4)]);
xlabel('\theta')
ylabel('Penalized Log-lkhd')
title('(d) L_2 Penalized Log-likelihood')

pause 
 
 
[a,in]=max(ell0);
ga0 = 0.01*1.05^in./stdx;

[a0,b0,c0,d0,e0,f0]=krigmle(x,y,ga0);

 

[a,in]=max(ell1);
ga1 = 0.01*1.05^in./stdx;

[lambda_sel1,rss1] = cvscad(x,y,ga1);

 
 

[a1,b1,c1,d1,e1,f1]=penkrigscad(x,y,lambda_sel1,ga1);

 
[a,in]=max(ell2);
ga2 = 0.01*1.05^in./stdx;

[lambda_sel2,rss2] = cvl1(x,y,ga2);

 

[a2,b2,c2,d2,e2,f2]=penkrigl1(x,y,lambda_sel2,ga2)

 

[a,in]=max(ell3);
ga3 = 0.01*1.05^in./stdx;

[lambda_sel3,rss3] = cvl2(x,y,ga3);

[a3,b3,c3,d3,e3,f3]=penkrigl2(x,y,lambda_sel3,ga3);

diary p12.out

[ga0,ga1,ga2,ga3]

[a0,a1,a2,a3]

[b0,b1,b2,b3]

[c0,c1,c2,c3]

[d0,d1,d2,d3]

[e0,e1,e2,e3]

[f0,f1,f2,f3]

[lambda_sel1,lambda_sel2,lambda_sel3]

rss1'
rss2'
rss3'

diary off



 

 