function data=unifgen;
u0=[
  4 28 13 25 30 16 12 12
 24 22 18 11 26 23  7  2
 14  7 11 28  7 24  6  6
 12 14 30  3  5 20  8 19
  3  9 19  5 20 28 24 13
  6  3 15 15 23 17  1 27
 10 17 25 27 25 26 17 22
 26 12  8  4 11 10  2 11
  7 25 24  7  3 12 14  5
 17 24  6 22 19 29  3 18
  8 13 21 21 12 14 30  1
 27 11 10 18 28 21 29 24
 21 20  7 26 22  6 23  4
 28 27 12  6 18  2 18 20
 18  8 29 14 29  8 20  7
 19  4  4  1 17 19 16  3
 16  2 14 12  1  4 28 15
 20 30 20 16  8  7  5 23
  9 15  1  8 27  5  9 17
 22  5 26  9 13 27 13 25
 11 29  3 13 14 25 27 10
 15 23 22  2 24 11 26 28
  5  6  5 23  6  9 22 21
  2 19 27 19 16  3  4  9
  1 21  9 10  9 22 19 26
 29  1 23 24 21 13 10 14
 13 10 17 30 15  1 15 30
 23 18  2 20  2 15 11 29
 30 16 16 17  4 30 21  8
 25 26 28 29 10 18 25 16];

u=(u0-0.5)/30;
u(:,1)=(0.15-0.05)*u(:,1)+0.05;
u(:,2)=(50000-100)*u(:,2)+100;
u(:,3)=(115600-63070)*u(:,3)+63070;
u(:,4)=(116-63.1)*u(:,4)+63.1;
u(:,5)=(1110-990)*u(:,5)+990;
u(:,6)=(820-700)*u(:,6)+700;
u(:,7)=(1680-1120)*u(:,7)+1120;
u(:,8)=(12045-9855)*u(:,1)+9855;

logrr=log(u(:,2)./u(:,1));
numer=2*pi*u(:,3).*(u(:,5)-u(:,6));
den=logrr.*(1+2*(u(:,7).*u(:,3))./(logrr.*(u(:,1).^2).*u(:,8))+u(:,3)./u(:,4));
yt=numer./den;
data=[[1:30]',u,yt];