load('piston12.m')
data = piston12;
y=data(:,8);
x=data(:,2:7);
 

[n,d]=size(x);

 
stdx=std(x)';
 

a=3.7;

ga1 = 0.01*1.05^45./stdx;

lam=[0.099,0.11,0.121]
aaa=[];
for l=1:3
lambda_sel1 = lam(l);

[mu,gamma0,beta0,si, ell,step]=penkrigscad(x,y,lambda_sel1,ga1);

aaa=[aaa,[mu;si;gamma0]];

load('piston100.m')
data = piston100;
z=data(:,2:7);
yp=data(:,8);

xt = x;

n=size(x,1);
q=2;

for k=1:100
   xv = z(k,:);

	R = exp(-sum( ((xt'-xv'*ones(1,n)).^q).*(gamma0*ones(1,n))));
   rss(k,l) = abs(yp(k) - mu - R*beta0);  
end;
end;
median(rss)

mean(rss)

figure(2)
 
subplot(221)
plot(rss(:,1),rss(:,2),'.',[0,4.5],[0,4.5],'-')
axis([0,4.5,0,4.5])
xlabel('AR of SCAD with \lambda_1')
ylabel('AR of SCAD with \lambda_2')
title('(a) Absolute Residual')

subplot(222)
plot(rss(:,1),rss(:,3),'.',[0,4.5],[0,4.5],'-')
axis([0,4.5,0,4.5])
xlabel('AR of SCAD with \lambda_1')
ylabel('AR of SCAD with \lambda_3')
title('(b) Absolute Residual')

subplot(223)
plot(rss(:,2),rss(:,3),'.',[0,4.5],[0,4.5],'-')
axis([0,4.5,0,4.5])
xlabel('AR of SCAD with \lambda_2')
ylabel('AR of SCAD with \lambda_3')
title('(c) Absolute Residual')



