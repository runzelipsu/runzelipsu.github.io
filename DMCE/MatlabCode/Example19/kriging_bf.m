function B = kriging_bf(x1,x2,theta)
[m1,n1] = size(x1);
[m2,n1] = size(x2);
% Calculate Correlation (Distance) Matrix
B = zeros(m1,m2);
for i=1:m1
    xt1 = x1(i,:);
    dx = (ones(m2,1)*theta).*(ones(m2,1)*xt1-x2).^2;
    if min(size(dx)) > 1
        d = sum(dx');
    else
        d = dx';
    end;
    B(i,:) = exp(-d);
end;

