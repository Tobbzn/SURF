% Saves Energy calculations to file.
function EnergySaver(freq_vec, Energy_vec, filename)
 figure
 title(filename)
 subplot(3,1,1)
 hold on
 plot(freq_vec,Energy_vec(:,end-2),'g')
 plot(freq_vec,Energy_vec(:,end-1),'r')
 plot(freq_vec,Energy_vec(:,end))
 plot(freq_vec,Energy_vec(:,end)+Energy_vec(:,end-2),'y')
 subplot(3,1,2)
 hold on
 plot(freq_vec,Energy_vec(:,end-2)./freq_vec,'g')
 plot(freq_vec,Energy_vec(:,end-1)./freq_vec,'r')
 plot(freq_vec,Energy_vec(:,end)./freq_vec)
   plot(freq_vec,Energy_vec(:,end-2)./freq_vec+Energy_vec(:,end)./freq_vec,'y')

 subplot(3,1,3)
 hold on
 plot(freq_vec,Energy_vec(:,end-2)./freq_vec.^2,'g')
 plot(freq_vec,Energy_vec(:,end-1)./freq_vec.^2,'r')
 plot(freq_vec,Energy_vec(:,end)./freq_vec.^2)
  plot(freq_vec,Energy_vec(:,end-2)./freq_vec.^2+Energy_vec(:,end)./freq_vec.^2,'y')
 xlabel('Frequency')
 ylabel('Energy')
% %savefig('test.fig')
save(['searchdata/' filename '.mat'], 'freq_vec','Energy_vec')
end