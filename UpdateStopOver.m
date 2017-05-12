function StopOver=UpdateStopOver(PopSize,Dimension)

[PopSize,Dimension]=size(PopSize);

for i=1:PopSize
        StopOver(i,:)=randperm(Dimension); %%Donor of searching for the stopover site
end
