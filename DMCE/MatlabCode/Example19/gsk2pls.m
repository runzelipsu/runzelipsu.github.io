function [pls,beta,s2,R, RInv] = gskpls(x,y,lambda,theta)
%function [y, Cx] = gskpls(datafile,lambda,theta)
[N,d] = size(x);
[like,beta,s2,R,RInv] = lhood(x,theta,y);
pls = like - N*sum(plambda(lambda,theta));

function p = penalty(lambda,theta)
p = sum(plambda(lambda,theta));

function p = plambda(lambda,theta)
p = lambda^2 - (abs(theta) - lambda).^2.*(abs(theta) < lambda);

function p = pdiv(lambda,theta)
a = 3.7;
eps = 1e-10;
p = lambda*(theta <= lambda) + 1/(a-1)*((a*lambda-theta).*(a*lambda > theta)).*(theta > lambda);

function p = pquad(lambda,theta,thetanew)
p = plambda(theta) + 0.5*(pdiv(lambda,theta)/abs(theta+eps)).*(thetanew.^2 - theta.^2);