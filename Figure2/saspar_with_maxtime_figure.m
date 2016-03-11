clf
load('saspar_with_maxtime_test1.mat')
plottable=[1,2,5,8];
marktable={'*-','o-','x-','*--','o--','+--','*-.','o-.'};
for i=plottable
plot(MaxTimeTable,record(i,:),marktable{i})
hold on
end
xlabel('max time allowed(s)')
ylabel('recovery rate (%)')
legend(num2str(xlengthtable(plottable)'),'Location','southeast')
title('Max Time Allowed vs Recovery Rate')
xlim([1,21])
ylim([10,110])