function Dij=Distance(xi,xj)

li=length(xi); 
lj=length(xj);
chromosome1=xi;
chromosome2=xj;
if rem(li,2)==0  %rem is function remainder after division of x by y
    M=(li-1)*li/2;
else
    M=(li-1)*(li+1)/2;
end
Dij=xj-xi;
s=0;
for i=1:li
    s=s+abs(Dij(i)); %Euclidean distance = ||xi-xj||
end
Dij=20*s/M;
