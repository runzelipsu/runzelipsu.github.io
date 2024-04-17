rand('seed',0);
randn('seed',0);

n=5000;
ni=50;
nj=80;
p=10;
nl=15;
displace=0.001;
nb=500;
ne=n;

xmu=0.0;
sigma=1;
sigma2=sigma*sigma;

du=1.0/(nj+1.0);
uj=[du:du:nj*du]';

beta00=sin(2*pi*uj);
beta11=4.0*(uj-0.5).*(uj-0.5);

dknot=1.0/(p+1.0);
uknot=[dknot:dknot:p*dknot]';

lambda1=logspace(log(0.001),log(1.5),nl)*log(nj)*sigma/sqrt(ni);
lambda2=logspace(log(0.001),log(2.5),nl)*log(nj)*sigma/sqrt(ni);

xp=[0:0.001:1]';
yp0=sin(2.0*pi*xp);
yp1=4*(xp-0.5).*(xp-0.5);

i=1;
beta0=beta00;
beta1=beta11;


gcvv0=[];
gcvv1=[];
gcv=[];
mse=[];
mmse=[];
varmse=[];
   
x=randn(ni,1);
eps=sigma*randn(ni,nj)+xmu;
y=ones(ni,1)*beta0';
y=y+x*beta1';
y=y+eps;

x=[ones(ni,1),x];
beta_hat=inv(x'*x)*x'*y;

mse(1,1)=(beta_hat(1,:)'-beta0)'*(beta_hat(1,:)'-beta0)/nj;
mse(1,2)=(beta_hat(2,:)'-beta1)'*(beta_hat(2,:)'-beta1)/nj;

%std(beta_hat'-[beta0,beta1])
beta0=beta_hat(1,:)';
beta1=beta_hat(2,:)';
   
%figure(1)
%plot(uj,beta0,'.',xp,yp0);
%figure(2)
%plot(uj,beta1,'.',xp,yp1);

for ii=1:nl
   %xlambda=lambda0(i);
   [betas, beta,gcv]=psr(uj,beta0,nj,uknot,p,lambda1(ii));
   gcvv0(ii)=gcv;
   [betas, beta,gcv]=psr(uj,beta1,nj,uknot,p,lambda2(ii));
   gcvv1(ii)=gcv; 
end;


[gcvmin0,indmin0]=min(gcvv0);
[beta0s, beta,gcv]=psr(uj,beta0,nj,uknot,p,lambda1(indmin0));
[gcvmin1,indmin1]=min(gcvv1);
[beta1s, beta,gcv]=psr(uj,beta1,nj,uknot,p,lambda2(indmin1));

%variance for smoothed coefficients
covt=sigma2*inv(x'*x);
for i=1:2
   varbi(i)=covt(i,i);
end;

[hg0]=psrvar(uj,beta0,nj,uknot,p,lambda1(indmin0));
[hg1]=psrvar(uj,beta1,nj,uknot,p,lambda1(indmin1));
covg0=varbi(1)*hg0*eye(nj)*hg0';
covg1=varbi(2)*hg1*eye(nj)*hg1';
for i=1,nj
   se0(i)=sqrt(covg0(i,i));
   se1(i)=sqrt(covg1(i,i));
end;
beta0s975=beta0s+1.96*se0;
beta0s025=beta0s-1.96*se0;
beta1s975=beta1s+1.96*se1;
beta1s025=beta1s-1.96*se1;

figure(1)
%plot(uj,beta0s,'k-', xp,yp0,'k:', uj,beta0s975,'k-.', uj,beta0s025,'k--');
%legend('Estimated', 'True', '97.5 percentile', '2.5 percentile',0);
subplot(2,1,1);
plot(uj,beta0s,'k-', xp,yp0,'k:');
%legend('Estimated', 'True',0);
axis([0 1 -1.5 1.5]);
xlabel('t');
ylabel('\beta_0');
title('(a) \beta_0');

%print -dpsc flmb0s1.eps;
   
   %figure(2)
   %plot(gcvv0,'*');
   %print -dpsc flmb0s4_gcv.eps;
%figure(3)
%plot(uj,beta1s,'k-', xp,yp1,'k:', uj,beta1s975,'k-.', uj,beta1s025,'k--');
%legend('Estimated', 'True', '97.5 percentile', '2.5 percentile',0);
subplot(2,1,2);
plot(uj,beta1s,'k-', xp,yp1,'k:');
%legend('Estimated', 'True',0);
axis([0 1 -0.5 1.5]);
xlabel('t');
ylabel('\beta_1');
title('(b) \beta_1');
%print -dpsc flmb1s1.eps;
print -dpsc flm-b-s1.eps;
   %figure(4)
   %plot(gcvv1,'*');
   %print -dpsc flmb1s4_gcv.eps;

%gcvmin0
%gcvmin1
   
mse(2,1)=(beta0s-beta00)'*(beta0s-beta00)/nj;
mse(2,2)=(beta1s-beta11)'*(beta1s-beta11)/nj;
   
for j=1:2
   for k=1:2
      varmse(j,k)=0.0;
      mmse(j,k)=0.0;
      mmse(j,k)=mmse(j,k)+(mse(j,k)-mmse(j,k))/i;
   end;
end;
mse11(i)=mse(1,1);
mse12(i)=mse(1,2);
mse21(i)=mse(2,1);
mse22(i)=mse(2,2);
mmse11(i)=mmse(1,1);
mmse12(i)=mmse(1,2);
mmse21(i)=mmse(2,1);
mmse22(i)=mmse(2,2);
varmse11(i)=varmse(1,1);
varmse12(i)=varmse(1,2);
varmse21(i)=varmse(2,1);
varmse22(i)=varmse(2,2);
rmse1(i)=mmse21(i)/mmse11(i);
rmse2(i)=mmse22(i)/mmse12(i);
rvarmse1(i)=varmse21(i)/varmse11(i);
rvarmse2(i)=varmse22(i)/varmse12(i);

   
for i=2:n
   beta0=beta00;
   beta1=beta11;

   gcvv0=[];
   gcvv1=[];
   gcv=[];
   
   x=randn(ni,1);
   eps=sigma*randn(ni,nj)+xmu;
   y=ones(ni,1)*beta0';
   y=y+x*beta1';
   y=y+eps;
   
   x=[ones(ni,1),x];
   beta_hat=inv(x'*x)*x'*y;
   
   mse(1,1)=(beta_hat(1,:)'-beta0)'*(beta_hat(1,:)'-beta0)/nj;
   mse(1,2)=(beta_hat(2,:)'-beta1)'*(beta_hat(2,:)'-beta1)/nj;
   
   %std(beta_hat'-[beta0,beta1]);
   beta0=beta_hat(1,:)';
   beta1=beta_hat(2,:)';
   
   %figure(1)
   %plot(uj,beta0,'.',xp,yp0);
   %figure(2)
   %plot(uj,beta1,'.',xp,yp1);
   
   for ii=1:nl
      %xlambda=lambda0(i);
      [betas, beta,gcv]=psr(uj,beta0,nj,uknot,p,lambda1(ii));
      gcvv0(ii)=gcv;
      [betas, beta,gcv]=psr(uj,beta1,nj,uknot,p,lambda2(ii));
      gcvv1(ii)=gcv; 
   end;
   
   [gcvmin0,indmin0]=min(gcvv0);
   [beta0s, beta,gcv]=psr(uj,beta0,nj,uknot,p,lambda1(indmin0));
   [gcvmin1,indmin1]=min(gcvv1);
   [beta1s, beta,gcv]=psr(uj,beta1,nj,uknot,p,lambda2(indmin1));
   %gcvmin0
   %gcvmin1
   
   mse(2,1)=(beta0s-beta00)'*(beta0s-beta00)/nj;
   mse(2,2)=(beta1s-beta11)'*(beta1s-beta11)/nj;
   
   for j=1:2
      for k=1:2
         varmse(j,k)=(1-1/(i-1))*varmse(j,k)+(mse(j,k)-mmse(j,k))^2/i;
         mmse(j,k)=mmse(j,k)+(mse(j,k)-mmse(j,k))/i;
         %         i,mmse(j,k),varmse(j,k)
%         mmse1(j,k,i)=mmse(j,k);
%         varmse1(j,k,i)=varmse(j,k);
      end;
   end;
   mse11(i)=mse(1,1);
   mse12(i)=mse(1,2);
   mse21(i)=mse(2,1);
   mse22(i)=mse(2,2);
   mmse11(i)=mmse(1,1);
   mmse12(i)=mmse(1,2);
   mmse21(i)=mmse(2,1);
   mmse22(i)=mmse(2,2);
   varmse11(i)=varmse(1,1);
   varmse12(i)=varmse(1,2);
   varmse21(i)=varmse(2,1);
   varmse22(i)=varmse(2,2);
	rmse1(i)=mmse21(i)/mmse11(i);
	rmse2(i)=mmse22(i)/mmse12(i);
	rvarmse1(i)=varmse21(i)/varmse11(i);
	rvarmse2(i)=varmse22(i)/varmse12(i);
   i
end;

save 'flm-mse' mse11 mse21 mse12 mse22;
save 'flm-mmse' mmse11 mmse21 mmse12 mmse22;
save 'flm-varmse' varmse11 varmse21 varmse12 varmse22;
save 'flm-ratio' rmse1 rmse2 rvarmse1 rvarmse2;

figure(4)
plot(mmse11(nb:ne));
ylabel('MMSE_0^1');
print -dpsc flm-mmse11.eps;

figure(5)
%plot(mmse(2,1,:));
plot(mmse21(nb:ne));
ylabel('MMSE_0^2');
print -dpsc flm-mmse21.eps;

figure(6)
plot(mmse12(nb:ne));
ylabel('MMSE_1^1');
print -dpsc flm-mmse12.eps;

figure(7)
plot(mmse22(nb:ne));
ylabel('MMSE_1^2');
print -dpsc flm-mmse22.eps;

figure(8)
plot(varmse11(nb:ne));
ylabel('VarMSE_0^1');
print -dpsc flm-varmse11.eps;

figure(9)
plot(varmse21(nb:ne));
ylabel('VarMSE_0^2');
print -dpsc flm-varmse21.eps;

figure(10)
plot(varmse12(nb:ne));
ylabel('VarMSE_1^1');
print -dpsc flm-varmse12.eps;

figure(11)
plot(varmse22(nb:ne));
ylabel('VarMSE_1^2');
print -dpsc flm-varmse22.eps;



figure(12)
plot(mse11(nb:ne));
ylabel('MSE_0^1');
print -dpsc flm-mse11.eps;

figure(13)
%plot(mmse(2,1,:));
plot(mse21(nb:ne));
ylabel('MSE_0^2');
print -dpsc flm-mse21.eps;

figure(14)
plot(mse12(nb:ne));
ylabel('MSE_1^1');
print -dpsc flm-mse12.eps;

figure(15)
plot(mse22(nb:ne));
ylabel('MSE_1^2');
print -dpsc flm-mse22.eps;


figure(16)
plot(rmse1(nb:ne));
ylabel('RMSE^1');
print -dpsc flm-rmse1.eps;

figure(17)
plot(rmse2(nb:ne));
ylabel('RMSE^2');
print -dpsc flm-rmse2.eps;

figure(18)
plot(rvarmse1(nb:ne));
ylabel('RVARMSE^1');
print -dpsc flm-rvarmse1.eps;

figure(19)
plot(rvarmse2(nb:ne));
ylabel('RVARMSE^2');
print -dpsc flm-rvarmse2.eps;

figure(20)
subplot(2,2,1);
plot([nb:1:ne],rmse1(nb:ne));
%title('a');
axis([nb ne  0.09 0.096]);
xlabel('k');
ylabel('MMSE^{(2)}_0/MMSE^{(1)}_0');
subplot(2,2,2);
plot([nb:1:ne],rmse2(nb:ne));
axis([nb ne  0.064 0.068]);
xlabel('k');
ylabel('MMSE^{(2)}_1/MMSE^{(1)}_1');
subplot(2,2,3);
plot([nb:1:ne],rvarmse1(nb:ne));
axis([nb ne  0.095 0.12]);
xlabel('k');
ylabel('VarMSE^{(2)}_0/VarMSE^{(1)}_0');
subplot(2,2,4);
plot([nb:1:ne],rvarmse2(nb:ne));
axis([nb ne  0.03 0.045]);
xlabel('k');
ylabel('VarMSE^{(2)}_1/VarMSE^{(1)}_1');
print -dpsc flm-ratio.eps;
