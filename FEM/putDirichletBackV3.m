function U = putDirichletBackV3(Uprev, lowerPlate, upperPlate, uzilow, uziup, x_plates, y_plates)

data = [uziup; uzilow; zeros(length(x_plates),1); zeros(length(y_plates),1)];
%[nodes, sortedindex] = sort([3*upperPlate; 3*lowerPlate;(3*x_plates-2); (3*y_plates-1)],'ascend');
%sorteddata = data(sortedindex);
U = zeros(length(Uprev)+length(data),1) + NaN;
U([3*upperPlate; 3*lowerPlate;(3*x_plates-2); (3*y_plates-1)]) = data;
%U(nodes+(0:(length(nodes)-1))') = sorteddata;
U(isnan(U)) = Uprev;

% for j = 1:length(nodes);
%     i = nodes(j);
%     temp = sorteddata(j);
%     N = length(U);
%     if i==1
%         U = [temp; U];
%     elseif i==N
%         U = [U; temp];
%     else
%         U = [U(1:(i-1)); temp; U(i:N)];
%     end
% end
 end