function Tp = ML_tleDistMatWin(A,n,p)
% ML_tleDistMatWin - averages temporally dinstance matrix A
% Written by Alan Peters with modification by Chad Jenkins and Michal
% Lewandowski

Tp=A;
[R C] = size(Tp);
D = zeros(R+n-1,C+n-1);
W = zeros(R+n-1,C+n-1);
Op = sparse(R,C);
% Tp = T;
% Op = O;
for k = 1:n
    for jj1=k:R+k-1
        for jj2=k:C+k-1
            D(jj1,jj2)=D(jj1,jj2)+Tp(jj1-k+1,jj2-k+1);
            if Op(jj1-k+1,jj2-k+1)==0
                W(jj1,jj2)=W(jj1,jj2)+1;
            end
        end
    end
%    D(k:R+k-1,k:C+k-1) = D(k:R+k-1,k:C+k-1) + Tp;
%    W(k:R+k-1,k:C+k-1) = W(k:R+k-1,k:C+k-1) + int32(not(Op));
    for jj1=1:R
        Tp(jj1,p-(k-1)) = 0;
    end
    for jj1=1:R
        Op(jj1,p-(k-1)) = 1;
    end
    for jj2=1:C
        Tp(p-(k-1),jj2) = 0;
    end
    for jj2=1:C
        Op(p-(k-1),jj2) = 1;
     end
end
    for jj1=1:R
        for jj2=1:C
            Tp(jj1,jj2)=D(jj1,jj2)/W(jj1,jj2);
        end
    end
%D = D(1:R,1:C);

end
