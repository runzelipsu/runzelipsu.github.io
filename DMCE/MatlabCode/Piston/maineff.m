load('piston12.m')
data = piston12;
y=data(:,8);
x=data(:,2:7);
[n,d]=size(x);

a =x(:,1);
[b,c]=sort(a);
d=y(c);

for k=1:6
   ay(k,1) = (d(2*k-1)+d(2*k))/2;
end;

subplot(3,2,1)
plot(a,y,'.',b(1:2:12),ay,'-')
xlabel('Clearance')
ylabel('Noise')

a =x(:,2);
[b,c]=sort(a);
d=y(c);

for k=1:6
   ay(k,1) = (d(2*k-1)+d(2*k))/2;
end;

subplot(3,2,2)
plot(a,y,'.',b(1:2:12),ay,'-')
xlabel('Press')
ylabel('Noise')


a =x(:,3);
[b,c]=sort(a);
d=y(c);

for k=1:6
   ay(k,1) = (d(2*k-1)+d(2*k))/2;
end;

subplot(3,2,3)
plot(a,y,'.',b(1:2:12),ay,'-')
xlabel('Skirt Length')
ylabel('Noise')


a =x(:,4);
[b,c]=sort(a);
d=y(c);

ay=zeros(3,1);
for k=0:2
   ay(k+1,1) = mean(d(4*k+1:4*(k+1)));
end;

subplot(3,2,4)
plot(a,y,'.',b(1:4:12),ay,'-')
xlabel('Skirt Profile')
ylabel('Noise')



a =x(:,5);
[b,c]=sort(a);
d=y(c);

ay=zeros(3,1);
for k=0:2
   ay(k+1,1) = mean(d(4*k+1:4*k+4));
end;

subplot(3,2,5)
plot(a,y,'.',b(1:4:12),ay,'-')
xlabel('Skirt Ovality')
ylabel('Noise')




a =x(:,6);
[b,c]=sort(a);
d=y(c);

ay=zeros(6,1);
for k=1:6
   ay(k,1) = (d(2*k-1)+d(2*k))/2;
end;

subplot(3,2,6)
plot(a,y,'.',b(1:2:12),ay,'-')

xlabel('Pin Offset')
ylabel('Noise')

