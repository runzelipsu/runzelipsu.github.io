function [lambda_sel,rss] = cvl1(x,y,ga0);
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

[n,d] = size(x);

q=2;
 
 

[n,d]=size(x);

la = 0.1*1.2.^[1:15]*sqrt(log(n)/n);
la = [0.03:0.05:0.75];

 
N=5;
for k=1:length(la)
   lambda0 = la(k);
 	rss(k)=0;	
    for i=1:N
        vd = [i:N:n];  % index for validation data set
        
       tr = [];
       for j=1:N  % To generate index for training data set
          if j~=i
             tr = [tr,[j:N:n]];
           end;
       end;
       xv = x(vd,:);  %  validation data set
       yv = y(vd);
       xt = x(tr,:); % Training data set
       yt = y(tr);
 
 		[mu_est,gamma_est, beta_est] = penkrigl1(xt,yt,lambda0,ga0);
 
		 % To compute R(x, x_i)	   
		 rss0 =[];
       for j=1:length(vd);
          xv0 = x(vd(j),:);
          yv0 = y(vd(j));
          n0=size(xt,1);
          R = exp(-sum( ((xt'-xv0'*ones(1,n0)).^q).*(gamma_est*ones(1,n0))));
          rss0 = [rss0,(yv0 - mu_est - R*beta_est)^2];
       end;
       
      rss(k) = rss(k) + sum(rss0);
   end;
end;
[cvs, index] = min(rss);
lambda_sel = la(index);

