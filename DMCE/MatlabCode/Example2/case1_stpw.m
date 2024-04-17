   a=[2   2   3   2   2   1   2   3   1.53   
   3   3   3   2   3   1   3   1   2.21 
   1   1   2   3   2   1   3   3   1.69 
   3   1   2   1   2   2   3   1   1.92 
   1   1   2   2   3   1   1   2   1.42 
   1   3   2   3   3   3   2   2   5.33 
   1   3   1   2   1   2   3   3   2.00 
   2   3   2   1   1   1   1   1   2.13 
   3   2   1   3   3   2   1   2   1.77 
   2   1   1   2   1   3   1   3   1.89 
   1   3   3   1   3   2   1   3   2.17 
   3   2   2   3   1   2   1   3   2.00 
   3   3   1   3   2   1   2   3   1.66 
   2   1   1   3   3   2   3   1   2.54 
   1   2   1   1   3   1   2   1   1.64 
   3   1   3   2   3   3   2   3   2.14 
   1   2   3   1   1   3   3   2   4.20 
   3   2   2   2   1   3   2   1   1.69 
   1   2   1   2   2   3   1   1   3.74 
   2   2   2   1   3   3   3   3   2.07 
   2   3   3   3   2   3   1   1   1.87 
   2   3   2   2   2   2   2   2   1.19 
   3   3   1   1   2   3   3   2   1.70 
   2   2   3   3   1   1   3   2   1.29 
   2   1   1   1   1   1   2   2   1.82 
   1   1   3   3   1   2   2   1   3.43 
   3   1   3   1   2   2   1   2   1.91 ];

x0=a(:,1:8);
y=a(:,9);

x=[x0,x0.^2];

for i=1:8
   for j=i+1:8
      x = [x,x0(:,i).*x0(:,j)];
   end;
end;
mx = mean(x);
sx = std(x);

[n,d] = size(x);

x = (x-ones(n,1)*mx)./(ones(n,1)*sx); %standardize x variables

y=y-mean(y);

r=log(n)/2;
bi=pinv(x'*x+eye(d,d)*r)*x'*y;
e=y-x*bi;

eff = trace(pinv(x'*x+eye(d,d)*r)*(x'*x))

s0 = sqrt(sum(e.^2)/(n-eff))

sigma=s0;

 
J1=stepwise(x,y,0.05,0.05);
x1=x(:,J1);
beta_s0 = (inv(x1'*x1)*x1'*y)';  % stepwise regression
d=size(x,2);
temp1 = zeros(1,d);
temp1(J1) = beta_s0;
beta_stpw = temp1

Jaic=stepwaic(x,y,sigma^2);
x1=x(:,Jaic);
beta_s0 = (inv(x1'*x1)*x1'*y)';  % stepwise regression
d=size(x,2);
temp1 = zeros(1,d);
temp1(Jaic) = beta_s0;
beta_aic = temp1;

Jbic=stepwbic(x,y,sigma^2);
x1=x(:,Jbic);
beta_s0 = (inv(x1'*x1)*x1'*y)';  % stepwise regression
temp1 = zeros(1,d);
temp1(Jbic) = beta_s0;
beta_bic = temp1;
   
Jric=stepwric(x,y,sigma^2);
x1=x(:,Jric);
beta_s0 = (inv(x1'*x1)*x1'*y)';  % stepwise regression
temp1 = zeros(1,d);
temp1(Jric) = beta_s0;
beta_ric = temp1;
   
Jphi=stepwphi(x,y,sigma^2);
x1=x(:,Jphi);
beta_s0 = (inv(x1'*x1)*x1'*y)';  % stepwise regression
temp1 = zeros(1,d);
temp1(Jphi) = beta_s0;
beta_phi = temp1;

[beta_stpw',beta_aic',beta_bic',beta_ric',beta_phi']
