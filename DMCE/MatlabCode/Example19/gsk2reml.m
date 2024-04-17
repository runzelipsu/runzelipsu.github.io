function [pls,beta,s2,R, RInv] = gsk2reml(x,y,lambda,theta)
%function [y, Cx] = gskpls(datafile,lambda,theta)
[N,d] = size(x);
[like,beta,s2,R,RInv] = lhoodreml(x,theta,y);
pls = like + 0.5*log(s2) -  0.5*log(ones(1,N)*RInv*ones(N,1));