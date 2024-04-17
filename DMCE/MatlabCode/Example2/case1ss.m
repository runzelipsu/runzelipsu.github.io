   a=[2   2   3   2   2   1   2   3   1.53   
   3   3   3   2   3   1   3   1   2.21 
   1   1   2   3   2   1   3   3   1.69 
   3   1   2   1   2   2   3   1   1.92 
   1   1   2   2   3   1   1   2   1.42 
   1   3   2   3   3   3   2   2   5.33 
   1   3   1   2   1   2   3   3   2.00 
   2   3   2   1   1   1   1   1   2.13 
   3   2   1   3   3   2   1   2   1.77 
   2   1   1   2   1   3   1   3   1.89 
   1   3   3   1   3   2   1   3   2.17 
   3   2   2   3   1   2   1   3   2.00 
   3   3   1   3   2   1   2   3   1.66 
   2   1   1   3   3   2   3   1   2.54 
   1   2   1   1   3   1   2   1   1.64 
   3   1   3   2   3   3   2   3   2.14 
   1   2   3   1   1   3   3   2   4.20 
   3   2   2   2   1   3   2   1   1.69 
   1   2   1   2   2   3   1   1   3.74 
   2   2   2   1   3   3   3   3   2.07 
   2   3   3   3   2   3   1   1   1.87 
   2   3   2   2   2   2   2   2   1.19 
   3   3   1   1   2   3   3   2   1.70 
   2   2   3   3   1   1   3   2   1.29 
   2   1   1   1   1   1   2   2   1.82 
   1   1   3   3   1   2   2   1   3.43 
   3   1   3   1   2   2   1   2   1.91 ];

x0=a(:,1:8);
y=a(:,9);

x=[x0,x0.^2];

for i=1:8
   for j=i+1:8
      x = [x,x0(:,i).*x0(:,j)];
   end;
end;
mx = mean(x);
sx = std(x);

[n,d] = size(x);

x = (x-ones(n,1)*mx)./(ones(n,1)*sx); %standardize x variables

ssto =y'*y;
 
sstoc = (y-mean(y))'*(y-mean(y));

ssr0 = n*mean(y)^2;

for k=1:d
   x1 = [ones(n,1),x(:,k)];
   ssr(k) = y'*x1*inv(x1'*x1)*x1'*y;
   sse(k) = ssto-ssr(k);
end;
ssr = ssr-ssr0;

for k=1:d
   if k~=6
      x1 = [ones(n,1),x(:,6),x(:,k)];
      ssr1(k) = y'*x1*inv(x1'*x1)*x1'*y;
      sse1(k) = ssto-ssr1(k);
   end;
end;
ssr1 = ssr1-ssr0-ssr(6);

for k=1:d
   if (k~=1)&(k~=6)
      x1 = [ones(n,1),x(:,1), x(:,6),x(:,k)];
      ssr2(k) = y'*x1*inv(x1'*x1)*x1'*y;
      sse2(k) = ssto-ssr2(k);
   end;
end;
ssr2 = ssr2-ssr0-ssr(6)-ssr1(1);

for k=8+1:8+8
      x1 = [ones(n,1),x(:,1), x(:,6),x(:,k)];
      ssr3(k) = y'*x1*inv(x1'*x1)*x1'*y;
      sse3(k) = ssto-ssr3(k);
 
end;
ssr3 = ssr3-ssr0-ssr(6)-ssr1(1);

for k=2*8+1:d
      x1 = [ones(n,1),x(:,1), x(:,6),x(:,k)];
      ssr4(k) = y'*x1*inv(x1'*x1)*x1'*y;
      sse4(k) = ssto-ssr4(k);
 
end;
ssr4 = ssr4-ssr0-ssr(6)-ssr1(1);

for k=1:d
   if (k~=1)&(k~=6)&(k~=21)
      x1 = [ones(n,1),x(:,1), x(:,6),x(:,21),x(:,k)];
      ssr5(k) = y'*x1*inv(x1'*x1)*x1'*y;
      sse5(k) = ssto-ssr5(k);
   end;
end;
ssr5 = ssr5-ssr0-ssr(6)-ssr1(1)-ssr4(21);

%%%%%%%%%%%%%%%%%%%%% Backward selection %%%%%%%%%%%%%%%%%

x2=x(:,[1,6,21]);


for k=1:3
      x1 = [ones(n,1),x2(:,k)];
      ssr6(k) = y'*x1*inv(x1'*x1)*x1'*y;
      sse6(k) = ssto-ssr6(k);
end;
ssr6 = ssr6-ssr0;


for k=1:3
   if k~=2
      x1 = [ones(n,1),x2(:,2),x2(:,k)];
      ssr7(k) = y'*x1*inv(x1'*x1)*x1'*y;
      sse7(k) = ssto-ssr7(k);
   end;
   
end;
ssr7 = ssr7-ssr0-ssr6(2);

for k=1:3
   if (k~=2)&(k~=3)
      x1 = [ones(n,1),x2(:,2),x2(:,3),x(:,k)];
      ssr8(k) = y'*x1*inv(x1'*x1)*x1'*y;
      sse8(k) = ssto-ssr8(k);
   end;
   
end;
ssr8 = ssr8-ssr0-ssr6(2)-ssr7(3);

