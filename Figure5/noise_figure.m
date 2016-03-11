clf
load noise.mat
tmprecord=[record(:,1:11,1);record(:,1:11,2)];
mark={'ro-','b*-','ks-','ro-.','b*-.','ks-.'};
for i=1:6
plot(snrtable(1:11),tmprecord(i,:)',mark{i})
hold on
end
legend('SASPAR: k=10','SASPAR: k=15','SASPAR: k=20','GESPAR: k=10','GESPAR: k=15','GESPAR: k=20','location','southeast')
xlabel('SNR (db)')
ylabel('Recovery Rate')
ylim([0,110])
title('Recover From Nosiy Measurements')