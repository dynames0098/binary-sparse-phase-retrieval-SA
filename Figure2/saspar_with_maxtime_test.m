
xlengthtable=[512:512:4096];
xsparsitytable=30:3:51;
supptable=xlengthtable/2;
MaxTimeTable=2:2:20;
repeat=100;
record=zeros(length(xlengthtable),length(MaxTimeTable));
for i3=1:repeat
    for i1=1:length(xlengthtable)
        for i2=1:length(MaxTimeTable)
            %%input
            xlength=xlengthtable(i1);
            xsparsity=xsparsitytable(i1);
            supp=supptable(i1);
            MaxTime=MaxTimeTable(i2);

            %generate sparse signal
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
            [~,loss]=sasper_maxtime(xlength,xsparsity,supp,y,G,Ginv,Gtrans,loss,MaxTime);
            if loss<1e-2
                record(i1,i2)=record(i1,i2)+1;
            end
            save('saspar_with_maxtime_test1.mat','record','xlengthtable','xsparsitytable','MaxTimeTable','repeat') 
            [i1,i2,i3]
        end
    end   
end
