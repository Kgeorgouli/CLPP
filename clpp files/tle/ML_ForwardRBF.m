function [P] = ML_ForwardRBF(Y,C,Bfor,sigma,rbfType)
%rbf - 0 thin plate, 1 - gaussian

    [P] = ML_ProjectRBF(Y,C,Bfor,sigma,rbfType);
    
end