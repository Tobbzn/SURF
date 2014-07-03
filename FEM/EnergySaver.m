% Saves Energy calculations to file.
function EnergySaver(freq_vec, Energy_vec, filename)
 figure
 title(filename)
 subplot(3,1,1)
 plot(freq_vec,Energy_vec(:,end))
 subplot(3,1,2)
 plot(freq_vec,Energy_vec(:,end)./freq_vec)
 subplot(3,1,3)
 plot(freq_vec,Energy_vec(:,end)./freq_vec.^2)
 xlabel('Frequency')
 ylabel('Energy')
% %savefig('test.fig')
save(['searchdata/' filename '.mat'], 'freq_vec','Energy_vec')
end