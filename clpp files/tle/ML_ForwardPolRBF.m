function [P] = ML_ForwardPolRBF(Y,Bfor,d)

    P=zeros(size(Y,1),d);
    for i=1:size(Y,1)
        phiX=Bfor*Y(i,:)';
        P(i,:)=phiX(end-d+1:end);
    end
end