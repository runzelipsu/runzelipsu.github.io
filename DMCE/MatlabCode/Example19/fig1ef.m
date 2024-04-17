function pls = fig1ab;

lambda=0; 
%x = [0 2 4 6 8 10]';
%x = (0:0.1:10)';
x=[0:0.5:10]'; 
y = sin(x); % - exp(x/100);
%y = sin(x);
 theta = (0.001:0.01:3)'; pmax = -10^6;

for i=1:length(theta);     
   ctheta = theta(i,:);     
   [p,beta,s2,R,RInv] = gsk2pls(x,y,lambda,ctheta);     
   pls(i) = p;
   if (p > pmax)         
      pmax = p;         
      bbeta = beta;         
      bs2 = s2;         
      btheta = ctheta;         
      bRInv = RInv;         
      %pause
     end;
end;
plsa = pls; 
[pmaxx,i] = max(plsa); 
subplot(2,2,3);
plot(theta,plsa);
theta(i);
title('(g) Log-Likelihood Function')
xlabel('\theta')
ylabel('Log-lkhd')
xp = [0:0.1:10]'; 
[d,d1] = size(x); 
for i=1:length(x)     
   x1 = x(i,:);     
   x2 = x;     
   r(i,:) = gsk_bf(x1,x2,btheta); 
end; 
gamma = regress(y-bbeta,r);

for i=1:length(xp);     
   x1 = xp(i,:);     
   x2 = x;     
   rc = gsk_bf(x1,x2,btheta);     
   yt(i) = bbeta +rc*gamma; 
end; 
%yp = sin(xp);  % - exp(xp/100) + 10;
yp = sin(xp); % - exp(xp/100);
subplot(2,2,4);
plot(x,y,'.',xp,yt,'-',xp,yp,'--'); 
title('(h) Prediction via Kriging')
xlabel('x')
ylabel('y')
axis([0,10,-1.1,1.1])
 