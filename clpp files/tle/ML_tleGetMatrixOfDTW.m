function [D,RD] = ML_tleGetMatrixOfDTW(X,tvwin,P)
% ML_tleGetMatrixOfDTW - calculates matrix of DTW distances
% Author   : Michal Lewandowski
%            Kingston University of London
%            Digital Imaging Research Centre
%            m.lewandowski@kingston.ac.uk

    N=size(X,1);
    D=zeros(N,N);
    RD=zeros(N,N);
    for i=1:N

        [iS,iE] = ML_tleGetWindowIndex(i,N,tvwin,P);
        [iRS,iRE] = ML_tleGetReverseWindowIndex(N-i+1,N,tvwin,P);
        for j=1:N

            [jS,jE] = ML_tleGetWindowIndex(j,N,tvwin,P);
            [jRS,jRE] = ML_tleGetReverseWindowIndex(N-j+1,N,tvwin,P);
            
            D(i,j) = ML_tleDTW(X(iS:iE,:),X(jS:jE,:));
            RD(N-i+1,N-j+1) = ML_tleDTW(X(iRS:iRE,:),X(jRS:jRE,:));
        end
%         disp(sprintf('DTW for %d ready',i)); Konstantia
    end
    D = D.*(1-eye(size(D)));
    RD = RD.*(1-eye(size(RD)));
end