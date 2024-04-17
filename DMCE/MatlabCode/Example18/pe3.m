x=[0:0.1:1]';
y=2*x.*cos(4*pi*x);

n=length(x);

xx=ones(n,1);

for k=1:9
       xx=[xx,x.^k];
end;

for k=1:10
for i=2:n-1
   x2=xx(i,1:k);
   y2=y(i);
%   if i==1
%      x1 = xx(2:n,1:k);
%      y1 = y(2:n);
%   end;
%   if i==n;
%      x1 = xx(1:n-1,1:k);
%      y1 = y(1:n-1);
%   end;
   
%   if (i~=1)&(i~=n)
      x1=[xx(1:i-1,1:k);xx(i+1:n,1:k)];
		y1=[y(1:i-1);y(i+1:n)];
%   end;
   
       b =inv(x1'*x1)*x1'*y1;
       pe00(i-1,k) = (y2-x2*b)^2;
 end;
 end;
 
 pe0 = mean(pe00)
 
a=[0:9];
figure(4)
 
plot(a,pe0,'-')

 
 