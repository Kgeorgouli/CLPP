function [D] = ML_tleGetMatrixOfDistMoments(DF,DB,tvwin,P,smoothSize)
% ML_tleGetMatrixOfDistMoments - combines matrices of DTW distances
% Author   : Michal Lewandowski
%            Kingston University of London
%            Digital Imaging Research Centre
%            m.lewandowski@kingston.ac.uk

    D=DF+DB;
    
    D(1:2*tvwin,1:P-2*tvwin) = DF(1:2*tvwin,1:P-2*tvwin).*2;
    D(1:P-2*tvwin,1:2*tvwin) = DF(1:P-2*tvwin,1:2*tvwin).*2;
    D(P-2*tvwin+1:P,2*tvwin+1:P) = DB(P-2*tvwin+1:P,2*tvwin+1:P).*2;
    D(2*tvwin+1:P,P-2*tvwin+1:P) = DB(2*tvwin+1:P,P-2*tvwin+1:P).*2;

    %clear DB DF DD1 DD2
 	D=ML_tleDistMatWin(D,smoothSize,P);
    D1=ML_tleDistMatWin(D(end:-1:1,end:-1:1),smoothSize,P);
    D=D+D1(end:-1:1,end:-1:1);
end