addpath c:\runze\rli\smooth
x=[0:0.1:1]';
x0=[0:0.01:1]';
y0=2*x0.*cos(4*pi*x0);
y=2*x.*cos(4*pi*x);

[yh0,grid]=gpnpr([x,y],0.04,[0,1,101],0,1);

[yh,grid1]=gpnpr([x,y],0.04,[0,1,11],0,1);


 

figure(1)

plot(x0,y0,'-',x,y,'o',grid,yh0,':',grid1,yh,'*')
legend('True', 'Sample','Fitted Curve','Fitted Values')
xlabel('x')
ylabel('y')

[yh2,grid2]=gpnpr([x,y],0.07,[0,1,101],0,1);

[yh3,grid3]=gpnpr([x,y],0.1,[0,1,101],0,1);

figure(2)

plot(x0,y0,'-',grid,yh0,'--',grid2,yh2,'-.',grid3,yh3,':')
legend('True', 'h=0.04','h=0.07','h=0.10')
xlabel('x')
ylabel('y')

