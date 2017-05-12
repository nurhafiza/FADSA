Gamma=50; %%light absorption (0-100)
p1=0.3*rand;
p2=0.3*rand;
MaxGen=10; 
PopSize=10; 
CompRun=10;

fireDiffDetail=fopen('fireDiffDetail','w');
fireDiffSim=fopen('fireDiffSim.txt','w');
fireDiffResults=fopen('fireDiffResults.txt','w');

ScoringMatrix=csvread('x60189_4.csv');
Dimension=size(ScoringMatrix,1); 

tic;
for i=1:CompRun
    fprintf('\nNo of Run %d\n',i);
    fprintf(fireDiffDetail,'\n====================');
    fprintf(fireDiffDetail,'\nNo of Run %d\n',i);
    fprintf(fireDiffDetail,'====================\n\n');
    fprintf(fireDiffSim,'\n====================');
    fprintf(fireDiffSim,'\nNo of Run %d\n',i);
    fprintf(fireDiffSim,'====================\n\n');

%% Initial population
Chrom=InitPop(PopSize,Dimension); 
ChromList=Chrom;

%% Initial solution with scoring value
%disp('An initial population of random solution: ');
InitialSolution=OutputSolution(Chrom(1,:));
InitialScoringValue=CalculateScore(ScoringMatrix,Chrom(1,:));

%% Optimization
ObjV=CalculateScore(ScoringMatrix,Chrom);
[preObjV,ObjVNo]=max(ObjV);

fprintf(fireDiffSim,'                     Best Solution\n\n');
for Gen=1:MaxGen
%    fprintf('\nIteration %d\n',Gen);
    fprintf(fireDiffDetail,'                    Gen      preObjV   CurrentBestSolution    TopCurrentBestSolution');
    %% Calculate fitness
%    ObjV=CalculateScore(ScoringMatrix,Chrom)
%    [preObjV,ObjVNo]=max(ObjV)
    
    for i=1:PopSize
        K=0;
        solution=zeros(PopSize,Dimension); 
        
        for j=1:i-1
            Dij=Distance(Chrom(i,:),Chrom(j,:)); 
            if Dij<=Gamma
                K=K+1; 
                solution(K,:)=Chrom(j,:); 
            end
        end
        
        for j=i+1:PopSize
            Dij=Distance(Chrom(i,:),Chrom(j,:)); 
             if Dij<=Gamma
                K=K+1;
                solution(K,:)=Chrom(j,:);
             end
        end
        
        if K==0
            solution1=zeros(1,Dimension);
        else
            solution1=ones(K,Dimension); 
        end
        
        for v=1:K
            solution1(v,:)=solution(v,:);
        end
        
        if K~=0
 %           fprintf('\nPopulation %d\n',i);
            fprintf(fireDiffDetail,'\nPopulation %d\n',i);
            ObjV1=CalculateScore(ScoringMatrix,solution1);
            [maxObjV1,ObjV1No]=max(ObjV1);
            maxChrom=solution1(ObjV1No,:);
%            ObjValue=CalculateScore(ScoringMatrix,Chrom(i,:))
            
            if preObjV<maxObjV1
%                fprintf('\nFA operation %d\n');
                fprintf(fireDiffDetail,'Firefly Algorithm');
                [Chrom(i,:),maxChrom]=ChangePlace(Chrom(i,:),maxChrom);
                CurrentBestSolution=maxObjV1;
                preObjV=CurrentBestSolution;           
            else 
                Donor=Chrom(randperm(PopSize),:);
                DonorScore=CalculateScore(ScoringMatrix,Donor);
                map=zeros(PopSize,Dimension);
                
                if rand<rand
                    if rand<p1
                        for m=1:PopSize
                            map(m,:)=rand(1,Dimension)<rand;
                        end
                    else
                        for m=1:PopSize
                            map(i,randi(Dimension))=1;
                        end
                    end
                else
                    for m=1:PopSize
                       map(m,randi(Dimension,1,ceil(p2*Dimension)))=1;
                    end 
                end
                
                Scale=4*randg;
                StopOver=PopSize+(Scale.*map).*(Donor-PopSize);
                StopOver=UpdateStopOver(StopOver,Dimension);
                StopOverObjV=CalculateScore(ScoringMatrix,StopOver);
                [MaxStopOverObjV,StopOverObjVNo]=max(StopOverObjV);
                
                if MaxStopOverObjV>preObjV 
                    fprintf(fireDiffDetail,'DS Algorithm     ');
%                    fprintf('\nDSA operation %d\n');
                    CurrentBestSolution=MaxStopOverObjV;
                    preObjV=CurrentBestSolution;
                else
                    fprintf(fireDiffDetail,'No Change        ');
%                    fprintf('\nNo Change %d\n');
                    CurrentBestSolution=preObjV;
                end
            end
            
                ShowNewObjValue(1,i)=CurrentBestSolution;                            
                TopCurrentBestSolution=max(ShowNewObjValue);
                fprintf(fireDiffDetail,'%5d  ---> %5d     ---> %5d         --->       %5d\n',Gen,preObjV,CurrentBestSolution,TopCurrentBestSolution);
        end                    
    end
                fprintf(fireDiffDetail,'--------------------------------------------------------------------------------------------------');
                BestSolution=max(ShowNewObjValue);
                globalmax(Gen)=BestSolution;
                iteration(Gen)=Gen;
                fprintf(fireDiffSim,'Iteration %5d   ---> %5d\n',Gen,BestSolution);
                clear ShowNewObjValue;
end
plot(iteration,globalmax)
xlabel('Iteration');
ylabel('Best Score');
hold on;

%ShowChromList=Chrom
%ScoreChromList=CalculateScore(ScoringMatrix,ShowChromList)
%OptimalSolution=OutputSolution(Chrom())
TopBestSolution=max(BestSolution);
fprintf(fireDiffSim,'\nTopBestSolution   ---> %5d\n',TopBestSolution);
fprintf(fireDiffResults,'%5d\n',TopBestSolution);
end
fclose(fireDiffDetail);
fclose(fireDiffSim);
fclose(fireDiffResults);
fireDiffResults=fopen('fireDiffResults.txt','r');
data=cell2mat(textscan(fireDiffResults,'%5d'));
data=dlmread('firediffresults.txt')
highestScore=max(data);
lowestScore=min(data);
avg=mean(data);
stdDev=std(data);
disp(['Best Scoring Value = ' num2str(highestScore)]);
disp(['Worst Scoring Value = ' num2str(lowestScore)]);
disp(['Average = ' num2str(avg)]);
disp(['Standard Deviation = ' num2str(stdDev)]);
fclose(fireDiffResults);
toc;