function [D] = ML_tleGetDTW_CTN(X,DX,tvwin,stdAcceptance,DTW1,DTW2,P,smoothSize,mintab)
% ML_tleGetDTW_CTN - constructs repetition heighbourhood graph
% Author   : Michal Lewandowski
%            Kingston University of London
%            Digital Imaging Research Centre
%            m.lewandowski@kingston.ac.uk

    N=size(X,1);
   

    if(isempty(mintab)==1 || ~exist('mintab','var'))
        [DF, DB] = ML_tleGetMatrixOfDistMoments1(DTW1,DTW2);
        DTW = ML_tleGetMatrixOfDistMoments(DF,DB,tvwin,P,smoothSize)
        [mintab] = ML_tleFindNeighbours(DTW,tvwin,P,stdAcceptance,smoothSize);
    end

    D=sparse(N,N);
    for i=1:N
        mt = mintab{i};   
        for m=1:size(mt,1)
%             [mS,mE] = ML_tleGetWindowIndex(mt(m,1),P,tvwin,P);
%             dist = DX(i,mS:mE)';
%             minD = min(dist);
%             ind = find(dist==minD);
%             if(length(ind)>1)
%             	ind = ind(1);
%             end
%             mt(m,2) = DX(i,mt(m,1));
%             mt(m,3) = mS+ind-1;
%             mt(m,4) = minD;
            
%             D(mt(m,1),i) = mt(m,2);
            D(mt(m,1),i) = DX(i,mt(m,1));
        end
        mintab{i} = mt;
    end          

    D(1:size(D,1)+1:end) = 1e-7;
    D = max(D, D');
    D = D .^ 2;
end