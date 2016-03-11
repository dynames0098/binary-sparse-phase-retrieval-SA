measurementstable=[500:50:1000];
xsparsitytable=[25:5:35];
repeat=100;
complexGaussian_result=zeros(length(xsparsitytable),length(measurementstable),2);
for i3=1:repeat
    for i1=1:length(xsparsitytable)
        for i2=1:length(measurementstable)
                
        %%input
        measurements=measurementstable(i2);
        xlength=measurementstable(i2);
        xsparsity=xsparsitytable(i1);
        supp=500;


        %generate sparse signal
        k=xsparsity;
        x=zeros(xlength,1);
        pp=randperm(supp-4);
        x(pp(1:xsparsity)+2)=1;


        %generate measurements
%         Gm=randn(measurements,xlength)+randn(measurements,xlength)*1j;
%         G=@(x)(Gm*x);
%         Ginv=@(x)(Gm\x);
%         Gtrans=@(x)(Gm'*x);
        G=@(x)fft(x);
        Ginv=@(x)ifft(x);
        Gtrans=@(x)ifft(x)*length(x);
        Gm=fft(eye(xlength));
        
         y=abs(G(x)).*abs(G(x));
         loss=@(x)(norm(abs(G(x)).^2-y)/norm(y));

        [~,loss]=sasper(xlength,xsparsity,supp,y,G,Ginv,Gtrans,loss);
        if loss<1e-2;
            complexGaussian_result(i1,i2,1)=complexGaussian_result(i1,i2,1)+1;
        end
        [~,l]=gesper(xlength,xsparsity,supp,G,Ginv,Gtrans,Gm,y);
        if l<1e-2;
            complexGaussian_result(i1,i2,2)=complexGaussian_result(i1,i2,2)+1;
        end
        progress=[i1,i2,i3]
        save('fewerMeasurements.mat','measurementstable','xsparsitytable','repeat','complexGaussian_result','progress')
        end
    end
end