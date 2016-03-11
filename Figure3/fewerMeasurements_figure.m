clf
load fewerMeasurements.mat

tmprecord=[complexGaussian_result(:,:,1);complexGaussian_result(:,:,2)];
mark={'ro-','b*-','ks-','ro-.','b*-.','ks-.'};
for i=1:6
plot(measurementstable(3:end),tmprecord(i,3:end)'/97*100,mark{i})
hold on
end
legend('SASPAR: k=25','SASPAR: k=30','SASPAR: k=35','GESPAR: k=25','GESPAR: k=30','GESPAR: k=35')
xlabel('number of measurements')
ylabel('Recovery Rate')
ylim([-10,110])
title('Recover From Fewer Measurements')