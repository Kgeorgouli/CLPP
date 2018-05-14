function [Dist]=ML_tleDTW(r,t)
% Continuous Dynamic Time Warping by Pau Mico 
%
% [Dist,D,k,w,rw,tw]=dtw(r,t,pflag)
%
% Dynamic Time Warping Algorithm
% Dist is unnormalized distance between t and r
% D is the accumulated distance matrix
% k is the normalizing factor
% w is the optimal path
% t is the vector you are testing against
% r is the vector you are testing
% rw is the warped r vector
% tw is the warped t vector
% pflag  plot flag: 1 (yes), 0(no)
%
% Version comments:
% rw, tw and pflag added by Pau Micó

[M,col]=size(r); %if (row > M) M=row; r=r'; end;
[N,col]=size(t); %if (row > N) N=row; t=t'; end;

d=ML_tleL2(r',t');

D=zeros(size(d));
D(1,1)=d(1,1);

for m=2:M
    D(m,1)=d(m,1)+D(m-1,1);
end
for n=2:N
    D(1,n)=d(1,n)+D(1,n-1);
end

for m=2:M
    for n=2:N
         D(m,n)=min(d(m,n)+D(m-1,n),min(d(m,n)+D(m-1,n-1),d(m,n)+D(m,n-1))); % this double MIn construction improves in 10-fold the Speed-up. Thanks Sven Mensing
    end
end
Dist=D(M,N)/(N+M);