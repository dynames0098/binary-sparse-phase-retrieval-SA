clf
load greedyVsSasPar2+
plot(xparsitytable,experiment_record1(2,:),'-*')
hold on
plot(xparsitytable,experiment_record1(1,:),'-o')
legend('SASPAR','Greedy')
title('SASPAR vs Greedy')
xlabel('sparsity k')
ylabel('recovery rate (%)')
xlim([22,50])
ylim([0,110])