rand('seed',1001);
randn('seed',1001);
 
ni=16;
nj=245;
p=10;
nl=15;
displace=0.001;

xmu=0.0;
sigma2=0.1;
sigma=sqrt(sigma2);

du=1.0/(nj+1.0);
uj=[du:du:nj*du]';

beta0=sin(2*pi*uj);
beta1=cos(2*pi*uj);
beta2=4.0*(uj-0.5).*(uj-0.5);
beta=[beta0,beta1,beta2];

dg=1.0/8.0;
gamma=ones(nj,1)*dg;
for i=2:7
gamma=[gamma,ones(nj,1)*dg*i];
end;

x=[randn(ni,2)<0.5];
x=[ones(ni,1),x];
z=+[randn(ni,7)<0.5];

eps=sigma*randn(nj,ni)+xmu;

y=beta*x'+gamma*z'+eps;

xt=[x,z];
betat=[beta,gamma];

betat_hat=y*xt*inv(xt'*xt);

std(betat_hat-betat)

%figure(1)
%plot(uj,betat(:,1:3),'.',uj,betat_hat(:,1:3));
%figure(2)
%plot(uj,betat(:,4:10),'.',uj,betat_hat(:,4:10));


dknot=1.0/(p+1.0);
uknot=[dknot:dknot:p*dknot]';

np=3;
for i=1:np
lambda(i,:)=logspace(log(0.001),log(1.5),nl)*log(nj)*std(betat_hat(:,i)-betat(:,i));
end;

gcvv0=[];
gcv=[];
bs=[];
%t1=(t-361)/amp;
xaxis=[];


for i=1:nl
xaxis(i)=i;
for j=1:np
[bs(:,j), beta,gcv]=csr1(uj,betat_hat(:,j),nj,uknot,p,lambda(j,i));
gcvv0(i,j)=gcv;
end;
end;

for j=1:np
[gcvmin0,indmin0]=min(gcvv0(:,j));
[bs(:,j), beta,gcv]=csr1(uj,betat_hat(:,j),nj,uknot,p,lambda(j,indmin0));
end;

%for i=4:10
%  i
%  sum(betat_hat(:,i))/nj
%  betat(1,i)
%  sum(betat_hat(:,i))/nj-betat(1,i)
%end;

%figure(1)
%plot(uj,bs(:,1:3),'-',uj,betat(:,1:3),':')

%figure(2)
%plot(xaxis,gcvv0(:,1),'-',xaxis,gcvv0(:,2:3),':')


%Now we got the initial value for gamma, begin to do the iteration
gamma0=betat_hat(:,4:10);
gamma0=repmat(mean(gamma0),[nj 1]);
gammao=gamma0;
gamma=gamma0
bso=bs;
error=1e-10;
maxerr=10.0;

while maxerr>error
  bso=bs;
  gamma=repmat(mean(gamma),[nj 1]);

  gamma0=gamma;
  
  ys=y-gamma0*z';
  beta_hat=ys*x*inv(x'*x);

  for i=1:np
    lambda(i,:)=logspace(log(0.001),log(1.5),nl)*log(nj)*std(beta_hat(:,i)-betat(:,i));
  end;

  gcvv0=[];
  gcv=[];
  bs=[];

  for i=1:nl
    for j=1:np
      [bs(:,j), beta,gcv]=csr1(uj,beta_hat(:,j),nj,uknot,p,lambda(j,i));
      gcvv0(i,j)=gcv;
    end;
  end;

  for j=1:np
    [gcvmin0,indmin0]=min(gcvv0(:,j));
    [bs(:,j), beta,gcv]=csr1(uj,beta_hat(:,j),nj,uknot,p,lambda(j,indmin0));
  end;

  ei=y-bs*x';

  gamma=ei*z*inv(z'*z);
  
%  maxerr=max(max(abs(gamma-gamma0)))
  maxerr=max(max(abs(bs-bso)))

end;

figure(1)
subplot(2,1,1);
plot(uj,beta0,'k:',uj,bs(:,1),'k-',uj,beta1,'k--',uj,bs(:,2),'k-.');
axis([0 1 -1 1.1]);
xlabel('t');
ylabel('\alpha_0, \alpha_1');
title('(a) \alpha_0 and \alpha_1');
%filename=strcat('e3alpha01.eps');
%print('-dpsc', filename);
subplot(2,1,2);
%figure(2)
plot(uj,beta2,'k:',uj,bs(:,3),'k-');
axis([0 1 -0.1 1.1]);
xlabel('t');
ylabel('\alpha_2');
title('(b) \alpha_2');
%filename=strcat('alpha2.eps');
%print('-dpsc', filename);
