function [ xbest,objectbest] = sasper_add_noise( xlength,xsparsity,supp,y,G,Ginv,Gtrans,loss,loss_accept)
%SASPERF Summary of this function goes here
%   Detailed explanation goes here
%%input
%parameters for sasper
%debug
countrestart=0;
MaxTime=20;
%debug
tic;
r = 0.99;
T = xlength/2;

k=xsparsity;
xlen=xlength;
setArea=supp;

%initialize
S=randperm(setArea);
S=S(1:k); %initialize S
xiter=zeros(xlen,1);
xiter(S)=1; %initialize x_iter
xbest=xiter;%initialize x_min
objectbest=loss(xiter); %initilaize F_min
objectlast=loss(xiter); %initilaize F_iter

while 1    
    T=r*T;
    Gx=G(xiter);
    grad1 = Gtrans(Gx.*(abs(Gx).*abs(Gx)-y)); 
    I1=grad1.*xiter.*(grad1>0); % compute I_1
    
    table=cumsum(I1);
    table=table/table(end);  %Distribution of i1
    i1=(find(rand<table,1,'first')); %select i1 randomly
    
    xMidIter=xiter;
    xMidIter(i1)=0; %hat x_iter
    
    Gx=G(xMidIter);
    grad2=Gtrans(Gx.*(abs(Gx).*abs(Gx)-y));  % nabla F(hat x_iter)
    I2=-grad2.*(1-xMidIter).*[ones(setArea,1);zeros(xlen-setArea,1)]; %I_2
    [~,index2]=max(I2);
    
    xnew=xMidIter;
    xnew(index2)=1;
    lossnew=loss(xnew);
    if lossnew<objectlast
       %newS=union(setdiff(S,index1),index2);
       %newSC=union(setdiff(nppc,index2),index1);
       xiter=xnew;
       objectlast=lossnew;
       if objectlast<objectbest
           xbest=xiter;
           objectbest=objectlast;
       end
    elseif rand(1)<exp((objectlast-lossnew/T))
        xiter=xnew;
        objectlast=lossnew;
    end
    if objectbest<loss_accept
        break;
    end
    if T<20
        %display('restart');
        T=xlength/2;
        S=randperm(setArea);
        S=S(1:k); %initialize S
        xiter=zeros(xlen,1);
        xiter(S)=1; %initialize x_iter
        countrestart=countrestart+1;
        objectlast=Inf;
        if toc>MaxTime
            display(['sasper Failed ',num2str(toc),' seconds passed'])
            break
        end
    end
end
end

