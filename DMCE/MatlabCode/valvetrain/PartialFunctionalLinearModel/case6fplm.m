clear
%rand('seed',1001)
%randn('seed',1001)

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
tt=a(nb:ne,1);
n=ne-nb+1;
nj=n;

y=a(nb:ne,2:p)'*100;

load('case6design.m')
x0=case6design;
% the 8th factor has 3 levels.
x01=[x0(:,2:7),(x0(:,8)==1),(x0(:,8)==2)];
[ni,np]=size(x01);

x01=[ones(ni,1),x01]; 
x=x01;
x(:,2)=x01(:,8);
x(:,3)=x01(:,9);
x(:,8)=x01(:,2);
x(:,9)=x01(:,3);

[ni,np]=size(x);

p=10;
nl=25;

betat_hat=inv(x'*x)*x'*y;
x00=x;

dknot=1.0/(p+1.0);
uknot=[dknot:dknot:p*dknot]';

np=3;
amp=360-90;
t1=(t-90)/amp;

z=x(:,4:9);
x=x(:,1:3);

%Initial value for gamma
gamma0=mean(betat_hat(4:9,:),2);
%Now we got the initial value for gamma, begin to do the iteration
gamma0
gamma=gamma0;
gammao=gamma0;
error=10^(-10);
maxerr=10.0;
ni=0;

while (maxerr>error)&(ni<200)
   gamma0=gamma;
   gamma0t=repmat(gamma0,[1 nj]);

   ys=y-z*gamma0t;
   beta_hat=inv(x'*x)*x'*ys;
   
   for i=1:np
      lambda(i,:)=logspace(log(0.0001),log(1.5),nl)*log(nj)*std(beta_hat(i,:));
   end;
   
   gcvv0=[];
   gcv=[];
   bs=[];
   bss=[];
   
   for i=1:nl
      for j=1:np
         beta_lv=beta_hat(j,:)';
         [bs, beta,gcv]=psr(t1,beta_lv,nj,uknot,p,lambda(j,i));
         gcvv0(i,j)=gcv;
      end;
   end;
   
   for j=1:np
      [gcvmin0,indmin0]=min(gcvv0(:,j));
      beta_lv=beta_hat(j,:)';
      [bs, beta,gcv]=psr(t1,beta_lv,nj,uknot,p,lambda(j,indmin0)/2);
%      [bs, beta,gcv]=psr(t1,beta_lv,nj,uknot,p,0.00000001);
      bss=[bss,bs];
   end;
   
   ei=y-x*bss';
   eia=mean(ei,2);
   
   gamma=inv(z'*z)*z'*eia;
   
   maxerr1=max(max(abs(gamma-gamma0)))
%   maxerr2=max(max(abs(bss-beta_hat')))
   
   maxerr=max(maxerr1,0)
   ni=ni+1
end;

figure(1)
subplot(2,2,1);
plot(tt,betat_hat(1,:),':',tt,bss(:,1),'-')
axis([90 360 -6 6]);
%legend('initial','smooth',0);
xlabel('Crank Angel');
ylabel('\beta_0');
title('(a) \beta_0');
%print -dpsc case6b0.eps;

subplot(2,2,2);
plot(tt,betat_hat(2,:),':',tt,bss(:,2),'-')
axis([90 360 -6 6]);
%legend('initial','smooth',0);
xlabel('Crank Angel');
ylabel('\beta_7');
title('(b) \beta_7');

%figure(2)
subplot(2,2,3);
plot(tt,betat_hat(3,:),':',tt,bss(:,3),'-')
axis([90 360 -6 6]);
%legend('initial','smooth',0);
xlabel('Crank Angel');
ylabel('\beta_8');
title('(c) \beta_8');
%print -dpsc case6b8.eps;


%figure(3)
%plot(tt,betat_hat(3,:),'k-',tt,bss(:,3),'k:')
%legend('initial','smooth',0);
%axis([90 360 -6 6]);
%xlabel('Crank Angel');
%print -dpsc case6b9.eps;

%figure(4)
%plot(gcvv0(:,1))
%figure(5)
%plot(gcvv0(:,2))
%figure(6)
%plot(gcvv0(:,3))

figure(7)
plot(tt,betat_hat(1:3,:),'k-',tt,betat_hat(4:9,:),'k:')
%legend('1st coefficient','8th coefficient','9th coefficient',0);
axis([90 360 -6 6]);
xlabel('Crank Angel');
ylabel('Initial LSE estimates for all coefficients');
%print -dpsc case6c_i.eps;


%%%% Calculate fitted values
b =gamma;
c=[bss,ones(n,1)*b'];
yhatss = c*x01';

y_flm=yhatss'/100;

 
yhatls = betat_hat'*x00';

y_ols=yhatls'/100;

y=y/100;

tt=tt';

%only plot run 1,5,9,13
figure(30)
subplot(2,2,1);
plot(tt,y(1,:),'k:',tt,y_flm(1,:),'k-');
axis([90 360 -0.06 .10]);
%    legend('CAE simulation','FLM estimate',0);
xlabel('Crank Angel');
ylabel('Motion Error');
title1=strcat('(a) Run 1')
    title(title1);
    subplot(2,2,2);
    plot(tt,y(5,:),'k:',tt,y_flm(5,:),'k-');
    axis([90 360 -.06 .10]);
%    legend('CAE simulation','FLM estimate',0);
    xlabel('Crank Angel');
    ylabel('Motion Error');
    title1=strcat('(b) Run 5')
    title(title1);
    subplot(2,2,3);
    plot(tt,y(9,:),'k:',tt,y_flm(9,:),'k-');
    axis([90 360 -.1 .10]);
%    legend('CAE simulation','FLM estimate',0);
    xlabel('Crank Angel');
    ylabel('Motion Error');
    title1=strcat('(c) Run 9')
    title(title1);
    subplot(2,2,4);
    plot(tt,y(13,:),'k:',tt,y_flm(13,:),'k-');
    axis([90 360 -.1 .10]);
%    legend('CAE simulation','FLM estimate',0);
    xlabel('Crank Angel');
    ylabel('Motion Error');
    title1=strcat('(d) Run 13')
    title(title1);

 
%only plot run 1,5,9,13 for OLS
figure(31)
    subplot(2,2,1);
    plot(tt,y(1,:),'k:',tt,y_ols(1,:),'k-');
    axis([90 360 -0.06 0.1]);
%    legend('CAE simulation','FLM estimate',0);
    xlabel('Crank Angel');
    ylabel('Motion Error');
    title1=strcat('(a) Run 1')
    title(title1);
    subplot(2,2,2);
    plot(tt,y(5,:),'k:',tt,y_ols(5,:),'k-');
    axis([90 360 -0.06 0.1]);
%    legend('CAE simulation','FLM estimate',0);
    xlabel('Crank Angel');
    ylabel('Motion Error');
    title1=strcat('(b) Run 5')
    title(title1);
    subplot(2,2,3);
    plot(tt,y(9,:),'k:',tt,y_ols(9,:),'k-');
    axis([90 360 -0.1 0.1]);
%    legend('CAE simulation','FLM estimate',0);
    xlabel('Crank Angel');
    ylabel('Motion Error');
    title1=strcat('(c) Run 9')
    title(title1);
    subplot(2,2,4);
    plot(tt,y(13,:),'k:',tt,y_ols(13,:),'k-');
    axis([90 360 -0.1 0.1]);
%    legend('CAE simulation','FLM estimate',0);
    xlabel('Crank Angel');
    ylabel('Motion Error');
    title1=strcat('(d) Run 13')
    title(title1);

