function [S,E] = ML_tleGetReverseWindowIndex(i,N,winSize,P)
% ML_tleGetReverseWindowIndex - takes index in relation to end of the sequence
% Author   : Michal Lewandowski
%            Kingston University of London
%            Digital Imaging Research Centre
%            m.lewandowski@kingston.ac.uk

    if(~exist('P','var'))
        P=N;
    end
    P=[1 P+1];
    [c,Pindex] = histc(1:1:N,[P]);
    
    %indeks z uwzglednieniem poczatku i konca
    E = i;
    if(E<P(Pindex(i))+2*winSize-1)
        E=P(Pindex(i))+2*winSize-1;
    end
    S = i-2*winSize+1;
    if(S<P(Pindex(i)))
        S = P(Pindex(i));
    end

end