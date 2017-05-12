function Chrom=InitPop(PopSize,Dimension)
 
Chrom=zeros(PopSize,Dimension); % For storing population
for i=1:PopSize
    Chrom(i,:)=randperm(Dimension); % Randomly generated initial population
end