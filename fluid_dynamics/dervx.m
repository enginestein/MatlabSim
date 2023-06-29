function D=dervx(N)

D(1,1)=1;
D(1,2)=2;

D(N,N)=1;
D(N,N-1)=-2;

 for i=2:N-1
        D(i,i-1)=-1;
        D(i,i)=0;
        D(i,i+1)=1;
end