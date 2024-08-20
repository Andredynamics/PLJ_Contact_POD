function [Phi,sigm,a]=POD(X,dt,r)

% sizeX=size(X);
% if sizeX(1)>sizeX(2)
%  
%     [Psi,lamb]= eig(X'*X);
%     sigm=sqrt(lamb);    
%     Phi=X*Psi/sigm;  
%     a=Phi'*X;
%     
% else

    [Phi,sigm,Psi] = svd(X,0);
    sigm=diag(sigm);
    a=Phi'*X;
% end

end