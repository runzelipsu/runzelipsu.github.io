function Cx = CorrMat(x,theta)
[m,n]=size(x);
Cx = ones(m,m);
for i=1:m-1
    x1 = x(i,:);
    x2 = x(i+1:m,:);
    d = gsk_bf(x1,x2,theta);
    Cx(i,i+1:m) = d;
    Cx(i+1:m,i) = d';
end;
    