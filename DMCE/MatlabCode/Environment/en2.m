a=[0.01   0.2  0.8 5.0 14.0  16.0    19.95   17.6  18.22  
0.05   2.0 10.0   0.1 8.0 12.0   22.09  22.85  22.62  
0.1 10.0   0.01  12.0   2.0 8.0  31.74  32.79  32.87  
0.2 18.0  1.0  0.8  0.4  4.0  39.37  40.65  37.87  
0.4  0.1 12.0  18.0  0.05   1.0  31.90  31.18  33.75  
0.8  1.0  0.05   4.0 18.0  0.4  31.14  30.66  31.18  
1.0 8.0  2.0  0.05  12.0   0.1  39.81  39.61  40.80  
2.0 16.0  14.0  10.0  5.0  0.01  42.48  41.86  43.79  
4.0  0.05   0.1  0.4  1.0 18.0    24.97  24.65  25.05 
5.0  0.8  4.0 16.0   0.2 14.0  50.29  51.22  50.54 
8.0 5.0 16.0   2.0  0.01  10.0  60.71  60.43  59.69 
10.0  14.0   0.2  0.01  16.0  5.0  67.01  71.99  67.12 
12.0   0.01  5.0 8.0 10.0   2.0  32.77  30.86  33.70  
14.0   0.4 18.0  0.2  4.0  0.8  29.94  28.68  30.66 
16.0   4.0  0.4 14.0   0.8  0.2  67.87  69.25  67.04 
18.0 12.0  8.0  1.0  0.1  0.05  55.56  55.28   56.52 
20.0 20.0  20.0 20.0  20.0 20.0  79.57 79.43 78.48];

y=mean(a(:,7:9)')';
x=a(:,1:6);
[n,d]=size(x);
mx=mean(x);
stdx=std(x);

x=(x-ones(n,1)*mx)./(ones(n,1)*stdx);
x=[ones(n,1),x];
b=inv(x'*x)*x'*y;
res = y-x*b;
rss2=sum(res.^2)

sigma = sqrt(sum(res.^2)/(n-7));
bse = sigma*sqrt(diag(inv(x'*x)));
[b,bse,2*(1-tcdf(abs(b./bse),n-7))]   

x1=x(:,1:3);


b1=inv(x1'*x1)*x1'*y;
res1 = y-x1*b1;
rss1=sum(res1.^2)

rss0=sum((y-mean(y)).^2) 

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

lambda_scad = gcv_scad([x,y],sigma^2);
[beta_scad,beta_scadstd] = scad(b,lambda_scad,[x,y],sigma^2);

 

lambda_lasso = gcv_lasso([x,y],sigma^2);
[beta_lasso,beta_lassostd] = lasso(b,lambda_lasso,[x,y],sigma^2);

lambda_scad
lambda_lasso

aa=[beta_stpw',beta_aic',beta_bic',beta_ric',beta_phi',beta_scad,beta_lasso]

mean((y*ones(1,7)-x*aa).^2)