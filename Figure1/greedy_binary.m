function [xbest,objectbest]=greedy_binary(xlength,xsparsity,supp,y,G,Ginv,Gtrans,loss)
tic
k=xsparsity;
xlen=xlength;
setArea=supp;
MaxTime=20;
%initialize
S=randperm(setArea);
S=S(1:k); %initialize S
xiter=zeros(xlen,1);
xiter(S)=1; %initialize x_iter
xbest=xiter;%initialize x_min
objectbest=loss(xiter); %initilaize F_min
objectlast=loss(xiter); %initilaize F_iter

while 1
    Gx=G(xiter);
    grad1 = Gtrans(Gx.*(abs(Gx).*abs(Gx)-y)); 
    I1=grad1.*xiter.*(grad1>0); % compute I_1
    [~,i1]=max(I1);
    xMidIter=xiter;
    xMidIter(i1)=0; %hat x_iter
    grad2=Gtrans(Gx.*(abs(G(xMidIter)).*abs(G(xMidIter))-y));  % nabla F(hat x_iter)
    I2=-grad2.*(1-xMidIter).*[ones(setArea,1);zeros(xlen-setArea,1)]; %I_2
    [~,index2]=max(I2);
    xnew=xMidIter;
    xnew(index2)=1;
    if loss(xnew)<objectlast
       xiter=xnew;
       objectlast=loss(xnew);
       if objectlast<objectbest
           xbest=xiter;
           objectbest=objectlast;
       end
    else
        S=randperm(setArea);
        S=S(1:k); %initialize S
        xiter=zeros(xlen,1);
        xiter(S)=1; %initialize x_iter
        objectlast=Inf;
       if toc>MaxTime
            display(['greedy binary Failed',num2str(toc),' seconds passed'])
            break
        end
    end
    if objectbest<1e-2
        break;
    end
end