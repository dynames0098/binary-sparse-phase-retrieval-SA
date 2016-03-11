function [x_new,F_new]=gespar_DGN(G,y)
%set parameters
Maxiter=100;
%inilization
A=real(G);
C=imag(G);
[~,n]=size(G);
x=randn(n,1);

F_iter=Inf;
x_iter=x;
step=0.5;
objectFun=@(x)(norm(abs(G*x).*abs(G*x)-y));
countIter=0;
while countIter<Maxiter
      step=min(1,2*step);
      countIter=countIter+1;
      B=(bsxfun(@times,A*x_iter,A)+bsxfun(@times,C*x_iter,C))*2;
      Gx=G*x_iter;
      b=y+abs(Gx).*abs(Gx);
      x_opt=B\b;
      x_new=x_iter+(x_opt-x_iter)*step;
      F_new=objectFun(x_new);
      while F_new > F_iter
            step=0.5*step;
            x_new=x_iter+(x_opt-x_iter)*step;
            F_new=objectFun(x_new);
      end
      if norm(x_new-x_iter)<1e-4
            return
      end
      x_iter=x_new;
      F_iter=F_new;
end
