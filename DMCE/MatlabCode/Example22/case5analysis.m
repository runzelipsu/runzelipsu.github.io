load('dmcase5.txt')
x=dmcase5;
x=[x(:,1:14),(x(:,15)==1),(x(:,15)==2),x(:,16:17)];
y=nvhresp;
[n,d]=size(x);
x=(x-ones(n,1)*mean(x))./(ones(n,1)*std(x));
ga0=0.01*ones(d,1).*(std(x)');
[mu_est,ga_est,beta_est,Psi2_est, ell,step]  = mvkrigmle(x,y,ga0)
