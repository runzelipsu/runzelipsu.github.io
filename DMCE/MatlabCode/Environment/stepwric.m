function J0 = stepwric(x,y,sigma2);
%
% stepwise variable selection procedure;
%
% Input: x  -- the design matrix
%        y  -- the response column vector
%        fin -- the significant level of F-enter;
%        fout -- the significant level of F-remove;
%
% Output: J0 = index of selected covariates.
%



[n,p]=size(x);
 
x0=[ ];
x1=x;
ric=[ ];
J0=[];
J1=[1:p];

% Forward selection to choose the 1st predictor

for j=1:p
   T=x(:,j);
   T0=y'*(eye(n)-T*inv(T'*T)*T')*y + 2*log(p)*size(T,2)*sigma2;
   ric=[ric;T0];
end;

[ric0,index]=min(ric);

J0 = [J0,index];

if index == 1
   J1=[2:p];
elseif index==p
   J1=[1:p-1];
else
   J1=[1:index-1,index+1:p];
end;

%%% End of first forward selection

%% Begin of stepwise selection %%

%% Forward selection

p0=length(J0);
p1=length(J1);



forward=1;
backward=1;

while ( p0<n ) & (length(p0)>0) & (length(p1)>0) & ((forward>0) | (backward>0)),
   
   
   x0 = x(:,J0);
   x1 = x(:,J1);
   p0=length(J0);
   p1=length(J1);

   ric = [ ];
   
   for j = 1:p1;
      T = [x0,x1(:,j)];
      T0 = y'*(eye(n)-T*inv(T'*T)*T')*y + 2*log(p) * size(T,2) * sigma2;
      ric=[ric;T0];
   end;
   
   [ric1,index] = min(ric);
   
    
     if ric1 < ric0
      J0=[J0,J1(index)];
      
      if index == 1;
         J1=J1(2:p1);
      elseif index==p
         J1=J1(1:p1-1);
      else
         J1=J1([1:index-1,index+1:p1]);
      end;
   else
      forward=0;
   end;
   
   x0=x(:,J0);
   x1=x(:,J1);
   p0=length(J0);
   p1=length(J1);
   ric0 = ric1;
   ric=[ ];
   
   for j=1:p0;
      if j==1
         T=x0(:,2:p0);
      elseif j==p;
         T=x0(:,1:p0-1);
      else
         T=x0(:,[1:j-1,j+1:p0]);
      end;
      
      T0=y'*(eye(n)-T*inv(T'*T)*T')*y + 2*log(p)* size(T,2) * sigma2;
      ric=[ric;T0];
   end;
   
   [ric1,index]=min(ric);
   
    
    
   if  ric1 < ric0
      J1=[J1,J0(index)];
      if index==1
         J0=J0(2:p0);
      elseif index==p0
         J0=J0(1:p0-1);
      else
         J0=J0([1:index-1,index+1:p0]);
      end;
   else
      backward=0;
   end;

   
end;

   
   
   
   
   
