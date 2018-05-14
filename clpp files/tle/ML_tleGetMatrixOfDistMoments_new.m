function [D] = ML_tleGetMatrixOfDistMoments(DD1,DD2,tvwin,P,smoothSize)
% ML_tleGetMatrixOfDistMoments - combines matrices of DTW distances
% Author   : Michal Lewandowski
%            Kingston University of London
%            Digital Imaging Research Centre
%            m.lewandowski@kingston.ac.uk

    if(size(DD1,1)~=size(DD1,2))
        DF=ML_L2(DD1',DD1',1);
        %DF=sparse(DF);
    else
        DF=DD1;
    end
    if(size(DD2,1)~=size(DD2,2))
        DB=ML_L2(DD2',DD2',1);  
        %DB=sparse(DB);
    else
        DB=DD2;
    end
    D=DF+DB;
    
    D(1:2*tvwin,1:P-2*tvwin) = DF(1:2*tvwin,1:P-2*tvwin).*2;
    D(1:P-2*tvwin,1:2*tvwin) = DF(1:P-2*tvwin,1:2*tvwin).*2;
    D(P-2*tvwin+1:P,2*tvwin+1:P) = DB(P-2*tvwin+1:P,2*tvwin+1:P).*2;
    D(2*tvwin+1:P,P-2*tvwin+1:P) = DB(2*tvwin+1:P,P-2*tvwin+1:P).*2;

    %clear DB DF DD1 DD2
 	DB=ML_tleDistMatWin(D,smoothSize,P);
    D1=DB(end:-1:1,end:-1:1);
    DF=ML_tleDistMatWin(D1,smoothSize,P);
    D=DB+DF(end:-1:1,end:-1:1);
end