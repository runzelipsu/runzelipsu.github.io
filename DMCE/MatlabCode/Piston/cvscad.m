function [lambda_sel,rss] = cvscad(x,y,ga0);
%
% lambda_sel = cvscad(x,y,vfold) is the leave-one-out cross-validation function
% for penalzed kriging method with SCAD penalty function
%
% Input: x -- x variables
%			y -- response variables
%			vfold -- v fold cross-validation
%
% Output: lambda_sel -- the selected tunging parameter
%

[n,d] = size(x)

q=2;
a=3.7;
 

[n,d]=size(x);

la = [0.03:0.05/2:0.75/2];

 
 

 
for k=1:length(la)
   lambda0 = la(k);
   rss(k)=0;
	 for i=1:n
       xv = x(i,:);  % To obtain the validation set
       yv = y(i);
      if i==1  % To obtain the training data set
         xt = x(2:n,:);
	      yt = y(2:n);
   	elseif i==n
      	xt = x(1:n-1,:);
	      yt = y(1:n-1);
   	else
      	xt = x([1:i-1,i+1:n],:);
	      yt = y([1:i-1,i+1:n]);
   	end;
      [mu_est,gamma_est, beta_est] = penkrigscad(xt,yt,lambda0,ga0);
      
      % To compute R(x,x_i)
      R = exp(-sum( ((xt'-xv'*ones(1,n-1)).^q).*(gamma_est*ones(1,n-1))));
      rss(k) = rss(k) + (yv - mu_est - R*beta_est)^2;
   end;
end;
[cvs, index] = min(rss);
lambda_sel = la(index);

      
      
      
   
      
   
