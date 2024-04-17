rand('seed',1001)
randn('seed',1001)

load('case6data.txt')
a=case6data;

clear case6data;
[n,p]=size(a);
t=a(:,1);
for i=1:n
  if t(i)==90
    nb=i
  elseif t(i)==360
    ne=i
  end
end;

t=a(nb:ne,1);
n=ne-nb+1;
nj=n;

y=a(nb:ne,2:p)*100;

load('case6design.m')
x0=case6design;
% the 8th factor has 3 levels.
x01=[x0(:,2:7),(x0(:,8)==1),(x0(:,8)==2)]; 
[ni,np]=size(x01);
 
x=[ones(ni,1),x01]; 
[ni,np]=size(x);

p=10;
nl=50;
displace=0.001;

betat_hat=y*x*inv(x'*x);

std(betat_hat)

dknot=1.0/(p+1.0);
uknot=[dknot:dknot:p*dknot]';

np=3;
for i=1:np
lambda(i,:)=logspace(log(0.0005),log(1.5),nl)*log(nj)*std(betat_hat(:,i));
end;

gcvv0=[];
gcv=[];
bs=[];
xaxis=[];
amp=360-90;
t1=(t-90)/amp;

for i=1:nl
xaxis(i)=i;
for j=1:np
[bs(:,j), beta,gcv]=csr1(t1,betat_hat(:,j),nj,uknot,p,lambda(j,i));
gcvv0(i,j)=gcv;
end;
end;

for j=1:np
[gcvmin0,indmin0]=min(gcvv0(:,j));
[bs(:,j), beta,gcv]=csr1(t1,betat_hat(:,j),nj,uknot,p,lambda(j,indmin0));
end;

%z=x(:,4:9);
%x=x(:,1:3);
z=x(:,2:7);
x = x(:,[1,8,9]);

 

 
%Now we got the initial value for gamma, begin to do the iteration
gamma0=sum(betat_hat(:,4:9))/nj;
gamma0
gamma=gamma0;
error=10^(-5);
maxerr=10.0;

while maxerr>error
  gamma0t=gamma0';
  for i=2:nj
    gamma0t=[gamma0t,gamma0'];
  end;
  gamma0t=gamma0t';
  ys=y-gamma0t*z';
  beta_hat=ys*x*inv(x'*x);

  for i=1:np
    lambda(i,:)=logspace(log(0.0001),log(1.5),nl)*log(nj)*std(beta_hat(:,i));
  end;

  gcvv0=[];
  gcv=[];
  bs=[];

  for i=1:nl
    for j=1:np
      [bs(:,j), beta,gcv]=csr1(t1,beta_hat(:,j),nj,uknot,p,lambda(j,i));
      gcvv0(i,j)=gcv;
    end;
  end;

  for j=1:np
    [gcvmin0,indmin0]=min(gcvv0(:,j));
    [bs(:,j), beta,gcv]=csr1(t1,beta_hat(:,j),nj,uknot,p,lambda(j,indmin0));
  end;

  ei=y-bs*x';

  eia=sum(ei)/nj;

  gamma=eia*z*inv(z'*z);

  maxerr=max(max(abs(gamma-gamma0)))
  gamma-gamma0

  gamma0=gamma;

end;
