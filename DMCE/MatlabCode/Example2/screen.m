

 
 
 
x=[x0,x0.^2];
x = (x-ones(n,1)*mean(x))./(ones(n,1)*std(x));  
d=size(x,2);
 
r=0;
bi=pinv(x'*x+eye(d,d)*r)*x'*y;
e=y-x*bi;

eff = trace(pinv(x'*x+eye(d,d)*r)*(x'*x))

s0 = sqrt(sum(e.^2)/(n-eff))

lambda= 1.175*sqrt(log(n))*s0/sqrt(n);

best= [];
bstd = [];
for k=12:12
 lambda= 4*k*0.025*sqrt(log(n))*s0/sqrt(n);
 lambda= 40*0.025*sqrt(log(n))*s0/sqrt(n);

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

tt = 4*[1:35]*sqrt(log(n))*s0/sqrt(n);

%plot(tt,gcv)
[bi(1:20),bscad0(1:20),std_scad0(1:20),bi(21:40),bscad0(21:40),std_scad0(21:40)]