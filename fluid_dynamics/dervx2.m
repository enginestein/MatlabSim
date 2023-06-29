function D2=dervx2(N)

D2(1,1)=2;
D2(1,2)=-5;
D2(1,3)=4;
D2(1,4)=-1;

D2(N,N)=2;
D2(N,N-1)=-5;
D2(N,N-2)=4;
D2(N,N-3)=-1;

 for i=2:N-1
        D2(i,i-1)=1;
        D2(i,i)=-2;
        D2(i,i+1)=1;
end