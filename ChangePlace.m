function [xi,xj]=ChangePlace(xi,xj)

[si,indexi]=sort(xi);
[sj,indexj]=sort(xj);
li=length(xi);
Tij=indexj-indexi;
r=rand(1,li);
alpha=0.50; % crossover happen
delta=0.95; % mutation happen
R=round(rand*2-1);
for i=1:li
    if r(i)<alpha
        indexi(i)=indexi(i);
    elseif r(i)<delta
        indexi(i)=indexj(i);
    else
        indexi(i)=indexj(i)+R;
    end
end
[si1,xi]=sort(indexi);
for i=1:li-1
    for j=i+1:li
        if si1(i)==si1(j)
            if Tij(i)<Tij(j)
                temp=xi(i);
                xi(i)=xi(j);
                xi(j)=temp;
            end
        end
    end
end