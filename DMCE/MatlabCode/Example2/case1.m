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

 

best= [];
bstd = [];
for k=1:25
 lambda= (1.05^(k-1))*0.2*sqrt(log(n))*s0/sqrt(n);
 %lambda= 40*0.025*sqrt(log(n))*s0/sqrt(n);

 [bscad0,std_scad0,gcv0,gcva0,gcvb0,effno0,effnoa0,effnob0]=scadls(x,y,lambda,bi);
 gcv(k,1)=gcv0;
 effno(k,1)=effno0;
 gcva(k,1)=gcva0;
 effnoa(k,1)=effnoa0;
 gcvb(k,1)=gcvb0;
 effnob(k,1)=effnob0;

 best = [best,bscad0];
 bstd =[bstd,std_scad0];
end;

tt = [1:100]*sqrt(log(n))*s0/sqrt(n);

plot(tt,gcv,'-',tt,gcva,':')
%[bi(1:20),bscad0(1:20),std_scad0(1:20),bi(21:40),bscad0(21:40),std_scad0(21:40)]



 