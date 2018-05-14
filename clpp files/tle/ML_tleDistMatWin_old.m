function D = ML_tleDistMatWin(A,n,p)
% ML_tleDistMatWin - averages temporally dinstance matrix A
% Written by Alan Peters with modification by Chad Jenkins and Michal
% Lewandowski

T=A;
[R C] = size(T);
D = zeros(R+n-1,C+n-1);
W = zeros(R+n-1,C+n-1);
O = ones(R,C);
Tp = T;
Op = O;
for k = 1:n
   D(k:R+k-1,k:C+k-1) = D(k:R+k-1,k:C+k-1) + Tp  ;
   W(k:R+k-1,k:C+k-1) = W(k:R+k-1,k:C+k-1) + Op ; 
   Tp(:,p-(k-1)) = 0;
   Tp(p-(k-1),:) = 0;
   Op(:,p-(k-1)) = 0;
   Op(p-(k-1),:) = 0;
end
D = D(1:R,1:C) ./ W(1:R,1:C);

end
