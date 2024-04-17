mu=[  56.72749999999997  56.25964746516307  56.51768178365280  56.53205821425562];

beta0= [
   0.02250000000003  -2.22970939596501  -0.55473229766242  -0.56483354136278
   0.92250000000006   3.46363103986515   1.12102122488036   1.80108665656768
  -2.75749999999997  -4.01560831896110  -1.49084773723854  -1.89518942412954
   2.04250000000004   1.34343449914094   2.43186542821453   2.38946688313278
  -0.38749999999997   0.72169111907329   0.57957638550628  -0.57195314265750
   0.12250000000003   0.00305383332517  -0.69641750016449  -1.26074270547951
  -0.04750000000025   1.19434791857423   0.11945126021751  -0.54206761348624
   1.72250000000003  -1.35638413163523   2.14639977109297   1.54297570863446
  -1.22749999999997  -0.51327274412463  -2.47237894227630  -1.74436596354077
  -3.95749999999997  -3.83668969824933  -3.26514718473346  -3.57609777565365
   0.63250000000003  -1.78982576957102   1.03186590671271   0.81279306947742
   2.91250000000003   7.01533164852743   3.91486827411173   3.60892784849752];

 
 
gamma0=[
     0.13971729636403   0.00082327821858   0.00166978761216   0.00377634396164
   1.63003512424697   0.00000018569523   0.00014179116251   0.02432604261031
   2.44505268637046   0.04268897371940   0.57785376347901   0.29086752208405
   4.09135569091184   0.00000056138195   0.00020223870912   0.03264154670866
   4.09135569091184   0.00000302741253   0.15014745443298   0.09798060789481
  12.22526343185230   4.62685914932326   0.01481465890209   0.25899413634225
];
 
 
load('piston12.m')
data = piston12;
y=data(:,8);
x=data(:,2:7);

xt = x;

minx = min(x);
maxx = max(x);
rx = maxx-minx;

K=10000;

z = rand(K,6);

n=size(x,1);
q=2;

 
f0=56.68436472184580;
D0= 2.83280675331207;
maineffdat;
%x0 = [0:0.01:1];

N=233;

x0 = mod([1:N]'*[1,89],N);
x0(N,:)=N;
x0 = (x0-0.5)/N;


N=length(x0);

dm=[1,2;
   1,3;
   1,4;
   1,5;
   1,6;
   2,3;
   2,4;
   2,5;
   2,6;
   3,4;
   3,5;
   3,6;
   4,5;
   4,6;
   5,6];
   
D2 = zeros(6,6);

for d=1:15
   d
   d1=dm(d,1);
   d2=dm(d,2);
      for l = 1:N
      	for k=1:K
        		xv = z(k,:).*rx+minx;
            xv(d1) = x0(l,1)*rx(d1)+minx(d1);
            xv(d2) = x0(l,2)*rx(d2)+minx(d2);  
      		R = exp(-sum( ((xt'-xv'*ones(1,n)).^q).*(gamma0(:,2)*ones(1,n))));
         	temp1(k) =  mu(2) + R*beta0(:,2);
         end;
         f1=interp1([0:0.01:1],f(:,d1),x0(l,1));
         f2=interp1([0:0.01:1],f(:,d2),x0(l,2));         
      	yhat2(l,1) = mean(temp1)-f0 -f1-f2;
      end;
      
      D2(d1,d2) = mean(yhat2.^2)*rx(d1)*rx(d2);
 
end;

diary globalsen2.out
format long;
D2
diary off
 


 