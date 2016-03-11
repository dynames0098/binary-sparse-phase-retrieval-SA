%GREEDY_binary vs saspar
repeat=100;
xparsitytable=[24:2:50];
experiment_record1=zeros(2,length(xparsitytable));
start_time=clock;
%%input
xlength=1024;
supp=512;
for i3=1:repeat
for i2=1:length(xparsitytable)
%generate sparse signal
xsparsity=xparsitytable(i2);
k=xsparsity;
x=zeros(xlength,1);
pp=randperm(supp-4);
x(pp(1:xsparsity)+2)=1;

%generate measurements
G=@(x)fft(x);
Ginv=@(x)ifft(x);
Gtrans=@(x)ifft(x)*length(x);
y=abs(G(x)).*abs(G(x));
loss=@(x)(norm(abs(G(x)).*abs(G(x))-y)/norm(y));

[~,l]=greedy_binary(xlength,xsparsity,supp,y,G,Ginv,Gtrans,loss);
if l<1e-4
    experiment_record1(1,i2)=experiment_record1(1,i2)+1;
end 
[~,l]=sasper(xlength,xsparsity,supp,y,G,Ginv,Gtrans,loss);
if l<1e-4
    experiment_record1(2,i2)=experiment_record1(2,i2)+1;
end 
[i2,i3]
end
save('greedyVsSasPar2+.mat','experiment_record1','xlength','supp','repeat','xparsitytable','start_time')
end
