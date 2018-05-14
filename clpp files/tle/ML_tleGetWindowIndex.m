function [S,E] = ML_tleGetWindowIndex(i,N,winSize,P)
% ML_tleGetWindowIndex - takes index in relation to start of the sequence
% Author   : Michal Lewandowski
%            Kingston University of London
%            Digital Imaging Research Centre
%            m.lewandowski@kingston.ac.uk

    if(~exist('P','var'))
        P=N;
    end
    [c,Pindex] = histc(1:1:N,[1,P+1]);

    S = i;
    if(S>P(Pindex(i))-2*winSize+1)
        S=P(Pindex(i))-2*winSize+1;
    end
    E = i+2*winSize-1;
    if(E>P(Pindex(i)))
        E = P(Pindex(i));
    end

end