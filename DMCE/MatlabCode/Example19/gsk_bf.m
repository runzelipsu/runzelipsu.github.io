function B = gsk_bf(x1,x2,theta)
[m1,n1] = size(x1);
[m2,n1] = size(x2);
% Calculate Correlation (Distance) Matrix
B = zeros(m1,m2);
for i=1:m1
    xt1 = x1(i,:);
    dx = repmat(theta,m2,1).*(repmat(xt1,m2,1)-x2).^2;
    d = sum(dx,2);
    B(i,:) = exp(-d');
end;

