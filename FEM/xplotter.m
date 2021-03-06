function [ xplot ] = xplotter( p,x,r,eps,n)

ballradius=max(p(:,3));
k=2*ballradius/n;
R = sqrt(p(:,1).^2+p(:,2).^2);
xplot = [];
for i=1:n
    K=ballradius-(i-1)*k;
    indices= find((abs(p(:,3)-K)<=k/2).*(R<r));
    if ~isempty(indices)
        xplot=[xplot;x(indices(1)*3,:)+K];
    end
end

    


end

