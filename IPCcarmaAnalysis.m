%% Carma video analysis
% Note: for trimming (i.e. trimvid) your n will shrink as you go above 600s
%% Properties that can be changed for trimming, time averaging, and fdr cutoff
time2condense=1; %Amount of time you want to average, make sure it is divisible by the trim time.
trimvid=900; %Specified time at which to trim for analysis.
cutoff=0.05;

load('carmaMeans.mat') %Might need to specify location if not in cd

%% Compile and average carma data
%To get an idea of the length of all the videos
% for i=1:61
%     timeLen(i,1) = length(carmaTimecourses{1,i});
% end
% timeLen=sort(timeLen);

timeCourses=nan(1163,61); %Needs to be the length of the longest video or analysis will be incorrect
for v=1:61
    lenVid=length(carmaTimecourses{1,v});
    timeCourses(1:lenVid,v)=carmaTimecourses{1:lenVid,v};
end
clear lenVid v

%To create a csv of just the averaged timecourses for each video
TCtable=writecarmacsv(timeCourses,controlExp);
writetable(TCtable,'carmaTC.csv') 

timeCourses=timeCourses(1:trimvid,:); %trim the timecourses  

%For comparing chunks of time
startCut=1;  
endCut=time2condense;
chunks=length(timeCourses)/time2condense;
for c=1:chunks
    for v=1:61
        timeChunk(c,v)=mean(timeCourses(startCut:endCut,v));
    end
    startCut=startCut+time2condense; 
    endCut=endCut+time2condense;
end
clear startCut endCut v c time2condense trimvid

controlTC=timeChunk(:,controlExp.experimental(:)==0);
expTC=timeChunk(:,controlExp.experimental(:)==1);

for ch=1:chunks
    [~,chunkPvals(ch,1),~,chunkStats(ch,:)] =ttest2(controlTC(ch,:),expTC(ch,:),'Vartype','unequal');  
end
clear ch

%% Multiple comparisons correction
% If we use Bonferroni here it is very strict 0.05/chunks=sig. cutoff
% bfCorr=0.05/chunks;
% sigMask=chunkPvals(:,1)<=bfCorr; %Nothing is sig.

% If we use FDR correction it seems to be the same
% [h, crit_p, adj_ci_cvrg, adj_p]=fdr_bh(chunkPvals,cutoff,'pdep','yes');
[h crit_p]=fdr_bky(chunkPvals,cutoff,'yes');
%% Plotting the significant data
for r=1:chunks
    meanExpTP(r,1)= nanmean(expTC(r,:));
    meanConTP(r,1)=nanmean(controlTC(r,:));
end
clear r
fdrMask=find(h==1); %Timepoints (i.e., rows that are sig.)

if ~isempty(fdrMask)
    conY=meanConTP(fdrMask); %Y coordinates for sig. timepoints 
    expY=meanExpTP(fdrMask); %Y coordinates for sig. timepoints 

    sigChunks=find(diff(fdrMask)>1); %are there gaps in sig? If so where
    numsigChunks=length(sigChunks)+1; %How many chunks of significance

    if numsigChunks > 1
        plot(meanConTP,'--', 'Color',[0, 0.4470, 0.7410])
        hold
        plot(meanExpTP,'--k')

        if numsigChunks == 2
            plot(fdrMask(1:sigChunks(1)),conY(1:sigChunks(1)),'-r')
            plot(fdrMask(sigChunks(1)+1:end),conY(sigChunks(1)+1:end),'-r')
            plot(fdrMask(1:sigChunks(1)),expY(1:sigChunks(1)),'-r')
            plot(fdrMask(sigChunks(1)+1:end),expY(sigChunks(1)+1:end),'-r')
        elseif numsigChunks == 3
            plot(fdrMask(1:sigChunks(1)),conY(1:sigChunks(1)),'-r')
            plot(fdrMask(sigChunks(1)+1:sigChunks(2)),conY(sigChunks(1)+1:sigChunks(2)),'-r')
            plot(fdrMask(sigChunks(2)+1:end),conY(sigChunks(2)+1:end),'-r')
            plot(fdrMask(1:sigChunks(1)),expY(1:sigChunks(1)),'-r')
            plot(fdrMask(sigChunks(1)+1:sigChunks(2)),expY(sigChunks(1)+1:sigChunks(2)),'-r')
            plot(fdrMask(sigChunks(2)+1:end),expY(sigChunks(2)+1:end),'-r')
        elseif numsigChunks == 4
            plot(fdrMask(1:sigChunks(1)),conY(1:sigChunks(1)),'-r')
            plot(fdrMask(sigChunks(1)+1:sigChunks(2)),conY(sigChunks(1)+1:sigChunks(2)),'-r')
            plot(fdrMask(sigChunks(2)+1:sigChunks(3)),conY(sigChunks(2)+1:sigChunks(3)),'-r')
            plot(fdrMask(sigChunks(3)+1:end),conY(sigChunks(3)+1:end),'-r')
            plot(fdrMask(1:sigChunks(1)),expY(1:sigChunks(1)),'-r')
            plot(fdrMask(sigChunks(1)+1:sigChunks(2)),expY(sigChunks(1)+1:sigChunks(2)),'-r')
            plot(fdrMask(sigChunks(2)+1:sigChunks(3)),expY(sigChunks(2)+1:sigChunks(3)),'-r')
            plot(fdrMask(sigChunks(1)+1:end),expY(sigChunks(3)+1:end),'-r')
        elseif numsigChunks == 5
            plot(fdrMask(1:sigChunks(1)),conY(1:sigChunks(1)),'-r')
            plot(fdrMask(sigChunks(1)+1:sigChunks(2)),conY(sigChunks(1)+1:sigChunks(2)),'-r')
            plot(fdrMask(sigChunks(2)+1:sigChunks(3)),conY(sigChunks(2)+1:sigChunks(3)),'-r')
            plot(fdrMask(sigChunks(3)+1:sigChunks(4)),conY(sigChunks(3)+1:sigChunks(4)),'-r')
            plot(fdrMask(sigChunks(4)+1:end),conY(sigChunks(4)+1:end),'-r')
            plot(fdrMask(1:sigChunks(1)),expY(1:sigChunks(1)),'-r')
            plot(fdrMask(sigChunks(1)+1:sigChunks(2)),expY(sigChunks(1)+1:sigChunks(2)),'-r')
            plot(fdrMask(sigChunks(2)+1:sigChunks(3)),expY(sigChunks(2)+1:sigChunks(3)),'-r')
            plot(fdrMask(sigChunks(3)+1:sigChunks(4)),expY(sigChunks(3)+1:sigChunks(4)),'-r')
            plot(fdrMask(sigChunks(4)+1:end),expY(sigChunks(4)+1:end),'-r')
        elseif numsigChunks == 6
            plot(fdrMask(1:sigChunks(1)),conY(1:sigChunks(1)),'-r')
            plot(fdrMask(sigChunks(1)+1:sigChunks(2)),conY(sigChunks(1)+1:sigChunks(2)),'-r')
            plot(fdrMask(sigChunks(2)+1:sigChunks(3)),conY(sigChunks(2)+1:sigChunks(3)),'-r')
            plot(fdrMask(sigChunks(3)+1:sigChunks(4)),conY(sigChunks(3)+1:sigChunks(4)),'-r')
            plot(fdrMask(sigChunks(4)+1:sigChunks(5)),conY(sigChunks(4)+1:sigChunks(5)),'-r')
            plot(fdrMask(sigChunks(5)+1:end),conY(sigChunks(5)+1:end),'-r')
            plot(fdrMask(1:sigChunks(1)),expY(1:sigChunks(1)),'-r')
            plot(fdrMask(sigChunks(1)+1:sigChunks(2)),expY(sigChunks(1)+1:sigChunks(2)),'-r')
            plot(fdrMask(sigChunks(2)+1:sigChunks(3)),expY(sigChunks(2)+1:sigChunks(3)),'-r')
            plot(fdrMask(sigChunks(3)+1:sigChunks(4)),expY(sigChunks(3)+1:sigChunks(4)),'-r')
            plot(fdrMask(sigChunks(4)+1:sigChunks(5)),expY(sigChunks(4)+1:sigChunks(5)),'-r')
            plot(fdrMask(sigChunks(5)+1:end),expY(sigChunks(5)+1:end),'-r')
        elseif numsigChunks == 7
            plot(fdrMask(1:sigChunks(1)),conY(1:sigChunks(1)),'-r')
            plot(fdrMask(sigChunks(1)+1:sigChunks(2)),conY(sigChunks(1)+1:sigChunks(2)),'-r')
            plot(fdrMask(sigChunks(2)+1:sigChunks(3)),conY(sigChunks(2)+1:sigChunks(3)),'-r')
            plot(fdrMask(sigChunks(3)+1:sigChunks(4)),conY(sigChunks(3)+1:sigChunks(4)),'-r')
            plot(fdrMask(sigChunks(4)+1:sigChunks(5)),conY(sigChunks(4)+1:sigChunks(5)),'-r')
            plot(fdrMask(sigChunks(5)+1:sigChunks(6)),conY(sigChunks(5)+1:sigChunks(6)),'-r')
            plot(fdrMask(sigChunks(6)+1:end),conY(sigChunks(6)+1:end),'-r')
            plot(fdrMask(1:sigChunks(1)),expY(1:sigChunks(1)),'-r')
            plot(fdrMask(sigChunks(1)+1:sigChunks(2)),expY(sigChunks(1)+1:sigChunks(2)),'-r')
            plot(fdrMask(sigChunks(2)+1:sigChunks(3)),expY(sigChunks(2)+1:sigChunks(3)),'-r')
            plot(fdrMask(sigChunks(3)+1:sigChunks(4)),expY(sigChunks(3)+1:sigChunks(4)),'-r')
            plot(fdrMask(sigChunks(4)+1:sigChunks(5)),expY(sigChunks(4)+1:sigChunks(5)),'-r')
            plot(fdrMask(sigChunks(5)+1:sigChunks(6)),expY(sigChunks(5)+1:sigChunks(6)),'-r')
            plot(fdrMask(sigChunks(6)+1:end),expY(sigChunks(6)+1:end),'-r')
        elseif numsigChunks == 8
            plot(fdrMask(1:sigChunks(1)),conY(1:sigChunks(1)),'-r')
            plot(fdrMask(sigChunks(1)+1:sigChunks(2)),conY(sigChunks(1)+1:sigChunks(2)),'-r')
            plot(fdrMask(sigChunks(2)+1:sigChunks(3)),conY(sigChunks(2)+1:sigChunks(3)),'-r')
            plot(fdrMask(sigChunks(3)+1:sigChunks(4)),conY(sigChunks(3)+1:sigChunks(4)),'-r')
            plot(fdrMask(sigChunks(4)+1:sigChunks(5)),conY(sigChunks(4)+1:sigChunks(5)),'-r')
            plot(fdrMask(sigChunks(5)+1:sigChunks(6)),conY(sigChunks(5)+1:sigChunks(6)),'-r')
            plot(fdrMask(sigChunks(6)+1:sigChunks(7)),conY(sigChunks(6)+1:sigChunks(7)),'-r')
            plot(fdrMask(sigChunks(7)+1:end),conY(sigChunks(7)+1:end),'-r')
            plot(fdrMask(1:sigChunks(1)),expY(1:sigChunks(1)),'-r')
            plot(fdrMask(sigChunks(1)+1:sigChunks(2)),expY(sigChunks(1)+1:sigChunks(2)),'-r')
            plot(fdrMask(sigChunks(2)+1:sigChunks(3)),expY(sigChunks(2)+1:sigChunks(3)),'-r')
            plot(fdrMask(sigChunks(3)+1:sigChunks(4)),expY(sigChunks(3)+1:sigChunks(4)),'-r')
            plot(fdrMask(sigChunks(4)+1:sigChunks(5)),expY(sigChunks(4)+1:sigChunks(5)),'-r')
            plot(fdrMask(sigChunks(5)+1:sigChunks(7)),expY(sigChunks(5)+1:sigChunks(6)),'-r')
            plot(fdrMask(sigChunks(6)+1:sigChunks(7)),expY(sigChunks(6)+1:sigChunks(7)),'-r')
            plot(fdrMask(sigChunks(7)+1:end),expY(sigChunks(7)+1:end),'-r')
        elseif numsigChunks == 9
            plot(fdrMask(1:sigChunks(1)),conY(1:sigChunks(1)),'-r')
            plot(fdrMask(sigChunks(1)+1:sigChunks(2)),conY(sigChunks(1)+1:sigChunks(2)),'-r')
            plot(fdrMask(sigChunks(2)+1:sigChunks(3)),conY(sigChunks(2)+1:sigChunks(3)),'-r')
            plot(fdrMask(sigChunks(3)+1:sigChunks(4)),conY(sigChunks(3)+1:sigChunks(4)),'-r')
            plot(fdrMask(sigChunks(4)+1:sigChunks(5)),conY(sigChunks(4)+1:sigChunks(5)),'-r')
            plot(fdrMask(sigChunks(5)+1:sigChunks(6)),conY(sigChunks(5)+1:sigChunks(6)),'-r')
            plot(fdrMask(sigChunks(6)+1:sigChunks(7)),conY(sigChunks(6)+1:sigChunks(7)),'-r')
            plot(fdrMask(sigChunks(7)+1:sigChunks(8)),conY(sigChunks(7)+1:sigChunks(8)),'-r')
            plot(fdrMask(sigChunks(8)+1:end),conY(sigChunks(8)+1:end),'-r')
            plot(fdrMask(1:sigChunks(1)),expY(1:sigChunks(1)),'-r')
            plot(fdrMask(sigChunks(1)+1:sigChunks(2)),expY(sigChunks(1)+1:sigChunks(2)),'-r')
            plot(fdrMask(sigChunks(2)+1:sigChunks(3)),expY(sigChunks(2)+1:sigChunks(3)),'-r')
            plot(fdrMask(sigChunks(3)+1:sigChunks(4)),expY(sigChunks(3)+1:sigChunks(4)),'-r')
            plot(fdrMask(sigChunks(4)+1:sigChunks(5)),expY(sigChunks(4)+1:sigChunks(5)),'-r')
            plot(fdrMask(sigChunks(5)+1:sigChunks(7)),expY(sigChunks(5)+1:sigChunks(6)),'-r')
            plot(fdrMask(sigChunks(6)+1:sigChunks(7)),expY(sigChunks(6)+1:sigChunks(7)),'-r')
            plot(fdrMask(sigChunks(7)+1:sigChunks(8)),expY(sigChunks(7)+1:sigChunks(8)),'-r')
            plot(fdrMask(sigChunks(8)+1:end),expY(sigChunks(8)+1:end),'-r')
        elseif numsigChunks == 10
            plot(fdrMask(1:sigChunks(1)),conY(1:sigChunks(1)),'-r')
            plot(fdrMask(sigChunks(1)+1:sigChunks(2)),conY(sigChunks(1)+1:sigChunks(2)),'-r')
            plot(fdrMask(sigChunks(2)+1:sigChunks(3)),conY(sigChunks(2)+1:sigChunks(3)),'-r')
            plot(fdrMask(sigChunks(3)+1:sigChunks(4)),conY(sigChunks(3)+1:sigChunks(4)),'-r')
            plot(fdrMask(sigChunks(4)+1:sigChunks(5)),conY(sigChunks(4)+1:sigChunks(5)),'-r')
            plot(fdrMask(sigChunks(5)+1:sigChunks(6)),conY(sigChunks(5)+1:sigChunks(6)),'-r')
            plot(fdrMask(sigChunks(6)+1:sigChunks(7)),conY(sigChunks(6)+1:sigChunks(7)),'-r')
            plot(fdrMask(sigChunks(7)+1:sigChunks(8)),conY(sigChunks(7)+1:sigChunks(8)),'-r')
            plot(fdrMask(sigChunks(8)+1:sigChunks(9)),conY(sigChunks(8)+1:sigChunks(9)),'-r')
            plot(fdrMask(sigChunks(9)+1:end),conY(sigChunks(9)+1:end),'-r')
            plot(fdrMask(1:sigChunks(1)),expY(1:sigChunks(1)),'-r')
            plot(fdrMask(sigChunks(1)+1:sigChunks(2)),expY(sigChunks(1)+1:sigChunks(2)),'-r')
            plot(fdrMask(sigChunks(2)+1:sigChunks(3)),expY(sigChunks(2)+1:sigChunks(3)),'-r')
            plot(fdrMask(sigChunks(3)+1:sigChunks(4)),expY(sigChunks(3)+1:sigChunks(4)),'-r')
            plot(fdrMask(sigChunks(4)+1:sigChunks(5)),expY(sigChunks(4)+1:sigChunks(5)),'-r')
            plot(fdrMask(sigChunks(5)+1:sigChunks(6)),expY(sigChunks(5)+1:sigChunks(6)),'-r')
            plot(fdrMask(sigChunks(6)+1:sigChunks(7)),expY(sigChunks(6)+1:sigChunks(7)),'-r')
            plot(fdrMask(sigChunks(7)+1:sigChunks(8)),expY(sigChunks(7)+1:sigChunks(8)),'-r')
            plot(fdrMask(sigChunks(8)+1:sigChunks(9)),expY(sigChunks(8)+1:sigChunks(9)),'-r')
            plot(fdrMask(sigChunks(9)+1:end),expY(sigChunks(9)+1:end),'-r')
        else
            disp('Sorry I have not made the loop for more than ten chunks of time')
        end
    elseif numsigChunks == 1
        plot(meanConTP,'--', 'Color',[0, 0.4470, 0.7410])
        hold
        plot(meanExpTP,'--k')
        plot(fdrMask,conY,'-r')
        plot(fdrMask,expY,'-r')        
    end
else
    disp('No significant chunks of time')
end

legend('Private','Public','Sig. Difference')
xlabel('Time (s)')
ylabel('Conflict (-100 to 100)')
%% Manually plot    
% plot(meanConTP,'--', 'Color',[0, 0.4470, 0.7410])
% hold
% plot(meanExpTP,'--k')
% plot(fdrMask,conY,'-r')
% plot(fdrMask,expY,'-r')
