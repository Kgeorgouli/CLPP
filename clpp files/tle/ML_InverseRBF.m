function [P] = ML_InverseRBF(X,C,Binv,sigma,rbfType)
%rbf - 0 thin plate, 1 - gaussian

    [P] = ML_ProjectRBF(X,C,Binv,sigma,rbfType);
    
end