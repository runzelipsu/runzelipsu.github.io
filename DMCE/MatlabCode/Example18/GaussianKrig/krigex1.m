 
lambda=0; 
x = [0:1:10]'/10;
y = 2*x.*cos(4*pi*x); % - exp(x/100);

theta = (0.5:0.05:10)'; 
pmax = -10^6;

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
subplot(3,2,1);
plot(theta,plsa);
title('(a) Log-Likelihood Function')
xlabel('\theta')
ylabel('Log-lkhd')
xp = (0:0.01:1)'; 
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
 
yp = 2*xp.*cos(4*pi*xp);  
subplot(3,2,2);
plot(xp,yp,'-',x,y,'o',xp,yt,':'); 
legend('True', 'Sample','Fitted Curve')
title('Prediction via Kriging')
xlabel('x')
ylabel('y')


figure(2)

plot(xp,yp,'-',x,y,'o',xp,yt,':'); 
legend('True', 'Sample','Fitted Curve')
title('Prediction via Kriging')
xlabel('x')
ylabel('y')

