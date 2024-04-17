load('boreholedata.dat')
data=boreholedata;
x0=data(:,2:9);
x0(:,1:2)=log(x0(:,1:2)) ;
y=data(:,10); % +0.01*randn(length(data),1);
y=log(y);
mx=mean(x0);
stdx=std(x0);
x=(x0-ones(length(y),1)*mx)./(ones(length(y),1)*stdx);

% To find quartiles

%sx=sort(x);
q1=prctile(x,15);
q2=median(x);
q3=prctile(x,85) ;

x2=[];
for i=1:size(x,2);
   for j=i:size(x,2)
      x2=[x2,x(:,i).*x(:,j)];  % interaction and quadratic terms
   end;   
end;

a1=x-ones(length(y),1)*q1;
a1=(a1.*(a1>0)).^2;

a2=x(:,1:2)-ones(length(y),1)*q2(1:2);
a2=(a2.*(a2>0)).^2;


a3=x-ones(length(y),1)*q3;
a3=(a3.*(a3>0)).^2;

d=[ones(length(y),1),x,x2,a1,a2,a3];

fin=0.1;
fout=0.1;


J0=stepwise(d,y,fin,fout)

d0=d(:,J0);

d1=d0;


%d1=d0(:,1:15);
beta0 = inv(d1'*d1)*d1'*y;

s0=sum((y-d1*beta0).^2)/(32-length(beta0))


lambda = gcv_scad([d0,y],s0)

betahat = inv(d0'*d0)*d0'*y;

[beta,std0] = scad(betahat,lambda,[d0,y],s0);

s0=sqrt(sum((y-d0*beta).^2)/(32-sum(beta~=0)))^2;

% Generate rand samples
N=10000;
u=rand(N,8);
u(:,1)=(0.15-0.05)*u(:,1)+0.05;
u(:,2)=(50000-100)*u(:,2)+100;
u(:,3)=(115600-63070)*u(:,3)+63070;
u(:,4)=(116-63.1)*u(:,4)+63.1;
u(:,5)=(1110-990)*u(:,5)+990;
u(:,6)=(820-700)*u(:,6)+700;
u(:,7)=(1680-1120)*u(:,7)+1120;
u(:,8)=(12045-9855)*u(:,1)+9855;

logrr=log(u(:,2)./u(:,1));
numer=2*pi*u(:,3).*(u(:,5)-u(:,6));
den=logrr.*(1+2*(u(:,7).*u(:,3))./(logrr.*(u(:,1).^2).*u(:,8))+u(:,3)./u(:,4));
yt=numer./den;

uu=u;
uu(:,1:2)=log(uu(:,1:2));
u0=(uu-ones(N,1)*mx)./(ones(N,1)*stdx);


u2=[];
for i=1:size(u0,2);
   for j=i:size(u0,2)
      u2=[u2,u0(:,i).*u0(:,j)];  % interaction and quadratic terms
   end;   
end;

us1=u0-ones(N,1)*q1;
us1=(us1.*(us1>0)).^2;

us2=u0(:,1:2)-ones(N,1)*q2(1:2);
us2=(us2.*(us2>0)).^2;


us3=u0-ones(N,1)*q3;
us3=(us3.*(us3>0)).^2;

ds=[ones(N,1),u0,u2,us1,us2,us3];

ds0=ds(:,J0);

yhat=exp(ds0*beta);

[J0',beta,std0]

mean((yt-yhat).^2)