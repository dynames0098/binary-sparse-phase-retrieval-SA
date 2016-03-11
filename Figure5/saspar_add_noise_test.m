xsparsitytable=[10,15,20];
snrtable=[4:1:20];
repeat=100;
record=zeros(length(xsparsitytable),length(snrtable),2);
start_time=clock;
for i3=1:repeat
    for i1=1:length(xsparsitytable)
        for i2=1:length(snrtable)
            %%input
            xlength=512;
            xsparsity=xsparsitytable(i1);
            supp=256;
            snrnow=snrtable(i2);
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

            loss=@(x)(norm(abs(G(x)).*abs(G(x))-y)/norm(y));
            lossaccept=norm(y_true-y)/norm(y)*1.1;
            trueloss=@(x)(norm(abs(G(x)).*abs(G(x))-y_true)/norm(y_true));

            [x,l]=sasper_add_noise(xlength,xsparsity,supp,y,G,Ginv,Gtrans,loss,lossaccept);
            % loss
            trueloss(x)
            if trueloss(x)<1e-4
                record(i1,i2,1)=record(i1,i2,1)+1;
            else 
                display('saspar failed')
            end
            
            [x,l]=gesper_add_noise(xlength,xsparsity,supp,G,Ginv,Gtrans,Gm,y,loss,lossaccept);
            x=round(x);
            % loss
            trueloss(x)
            if trueloss(x)<1e-4
                record(i1,i2,2)=record(i1,i2,2)+1;
           else 
                display('gespar failed')
            end
            progress=[i1,i2,i3]
            
        end
    end
    end_time=clock;
    save('noise.mat','record','xsparsitytable','snrtable','repeat','start_time','end_time','progress')
end
