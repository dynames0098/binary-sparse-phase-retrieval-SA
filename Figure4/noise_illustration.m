clf

xlength=512;
xsparsity=10;
supp=256;
snrnow=6;
%generate sparse signal
k=xsparsity;
x=zeros(xlength,1);
pp=randperm(supp-4);
x(pp(1:xsparsity)+2)=1;

%generate measurements
Gm=fft(eye(xlength));
G=@(x)fft(x);
Ginv=@(x)ifft(x);
Gtrans=@(x)ifft(x)*length(x);
y=abs(G(x)).*abs(G(x));
y_true=y;
%add noise
y=awgn(y_true,snrnow,'measured');
plot(y_true,'o-')
hold on
plot([y],'*--')


plot([100,200],[0,0],'--')
ylim([-20,60])
xlim([100,200])
legend('accurate measurements','noisy measurements (4db)','location','southeast')
title('Noisy Measurements')