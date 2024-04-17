clear
 

load('case6data.txt')
a=case6data;

clear case6data;
[n,p]=size(a);
t=a(:,1);
for i=1:n
  if t(i)==360
    nb=i
  elseif t(i)==450
    ne=i
  end
end;

t=a(nb:ne,1);
tt=a(nb:ne,1);
n=ne-nb+1;
nj=n;

y=a(nb:ne,2:p)'*100;

load('case6design.m')
x0=case6design;
% the 8th factor has 3 levels.
x01=[x0(:,2:7),(x0(:,8)==1),(x0(:,8)==2)]; 

xk = x01;

[ni,np]=size(x01);

x01=[ones(ni,1),x01]; 
x=x01;

[ni,np]=size(x);

p=10;
nl=25;

h1=inv(x'*x)*x';
betat_hat=h1*y;

dknot=1.0/(p+1.0);
uknot=[dknot:dknot:p*dknot]';

amp=450-360;
t1=(t-360)/amp;

for i=1:np
   lambda(i,:)=logspace(log(0.0001),log(1.5),nl)*log(nj)*std(betat_hat(i,:));
end;
   
gcvv0=[];
gcv=[];
bs=[];
bss=[];
bsst=[];
bssb=[];
   
for i=1:nl
   for j=1:np
      beta_lv=betat_hat(j,:)';
      [bs, beta,gcv]=psr(t1,beta_lv,nj,uknot,p,lambda(j,i));
      gcvv0(i,j)=gcv;
   end;
end;

%calculate the confidence interval for the estimated coefficients bss
yhatt=x*betat_hat;
hsigma2=0;
for i=1:np
   hsigma2=hsigma2+(y(i,:)-yhatt(i,:))*(y(i,:)-yhatt(i,:))';
end;
hsigma2=hsigma2/(ni-np)/(nj-1);
covt=hsigma2*inv(x'*x);
for i=1:np
   varbi(i)=covt(i,i);
end;

for j=1:np
   [gcvmin0,indmin0]=min(gcvv0(:,j));
   beta_lv=betat_hat(j,:)';
   [bs, beta,gcv]=psr(t1,beta_lv,nj,uknot,p,lambda(j,indmin0));
   [hg]=psrvar(t1,beta_lv,nj,uknot,p,0.00000001);
   covg=varbi(j)*hg*eye(nj)*hg';
   bst=bs+1.96*sqrt(diag(covg));
   bsb=bs-1.96*sqrt(diag(covg));
   bss=[bss,bs];
   bsst=[bsst,bst];
   bssb=[bssb,bsb];
end;


%calculate the Functional Linear Model estimate using X\hat{\beta}
y_flm=x*bss';
%y=y/100;
%y_flm=y_flm/100;
y_ols=x*betat_hat;
%y_ols=y_ols/100;

yrflm = y-y_flm;
yrols = y-y_ols;

muf=[];
sigma2f=[];
gaf=[];

mul=[];
sigma2l=[];
gal=[];


for i=1:size(yrflm,2)
  
   [mu_est,sigma2_est,ga_est,step] = krig(xk,yrflm(:,i));
   muf=[muf,mu_est];
   sigma2f = [sigma2f,sigma2_est];
   gaf = [gaf,ga_est];
end;

 

fid=fopen('gaflmkirg.m','w');
fprintf(fid,'%7.4f %7.4f %7.4f %7.4f %7.4f %7.4f %7.4f %7.4f\n', gaf);
fclose(fid)

 

fid=fopen('sigma2flmkirg.m','w');
fprintf(fid,'%7.4f\n', sigma2f);
fclose(fid)

 

fid=fopen('muflmkirg.m','w');
fprintf(fid,'%7.4f\n', muf);
fclose(fid)

 



