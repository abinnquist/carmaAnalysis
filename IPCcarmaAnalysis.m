%% Carma video analysis
% Note: for trimming (i.e. trimvid) your n will shrink as you go above 600s
carmaPath=uigetdir('','Choose Data Directory');

%% Properties that can be changed for trimming, time averaging, and fdr cutoff
numConvos=61;
trimvid=900; %Specified time at which to trim for analysis.
cutoff=0.05;
time2condense=10; %Amount of time you want to average, make sure it's divisible by trim time.

overrideH = 0.0;
conflictSort=4;

%% Analyses to perform
condense=0;         % If you want to average chunks of time instead of each time point
halves=1;           % Compare the first half of the conversation to the second
medianSplit=0;

%% Load in external data 
load(strcat(carmaPath,filesep,'carmaMeans.mat'))
IPCdata=readtable(strcat(carmaPath,filesep,'IPCdata_complete.csv'));

%% Compile and average carma data
timeCourses=nan(1163,61); %Needs to be the length of the longest video or analysis will be incorrect
for v=1:61
    lenVid=length(carmaTimecourses{1,v});
    timeCourses(1:lenVid,v)=carmaTimecourses{1:lenVid,v};
end

%trim the timecourses  
timeChunk=timeCourses(1:trimvid,:); 

if condense
    %For comparing chunks of time
    startCut=1;  
    endCut=time2condense;
    chunks=length(timeCourses)/time2condense;
    for c=1:chunks
        for v=1:numConvos
            timeChunk(c,v)=mean(timeCourses(startCut:endCut,v));
        end
        startCut=startCut+time2condense; 
        endCut=endCut+time2condense;
    end
    clear startCut endCut v c time2condense
end

controlTC=timeChunk(:,controlExp.experimental(:)==0);
expTC=timeChunk(:,controlExp.experimental(:)==1);

%% t-test between first half and last half of the conversation
% Meaning into a single timecourse first, then doing analysis
if halves
    avgTimecourse = mean(timeChunk,2,'omitnan');
    firsthalf = avgTimecourse(1:trimvid/2);
    lasthalf = avgTimecourse((trimvid/2)+1:end);

    f = categorical({'First Half','Second Half'}); % Predefine this for easier plotting
    f = reordercats(f,{'First Half','Second Half'});
    [~, p, ~, stats] = ttest2(firsthalf, lasthalf)
    figure()
    [t, mean_a, mean_b, se_a, se_b] = ttest2_barplot(firsthalf,lasthalf,f);
    ylabel('Conflict')

    % Using all datapoints for all timecourses
    firstAll = reshape(timeChunk(1:trimvid/2,:),1,[]);
    lastAll = reshape(timeChunk((trimvid/2)+1:end,:),1,[]);

    [~, p, ~, stats] = ttest2(firstAll, lastAll)
    figure()
    [t, mean_a, mean_b, se_a, se_b] = ttest2_barplot(firstAll,lastAll,f);
    ylabel('Conflict')
end

%% Median split of first half last half (2-way ANOVA)
if medianSplit
    first450SplitCon = mean(controlTC(1:450,:));
    first450SplitExp = mean(expTC(1:450,:));
    last450SplitCon = mean(controlTC(451:end,:));
    last450SplitExp = mean(expTC(451:end,:));
    anovaSplit = [first450SplitCon first450SplitExp last450SplitCon last450SplitExp];
    groupFirst = [repmat({'First'},1,length(first450SplitCon)) repmat({'First'},1,length(first450SplitExp))...
         repmat({'Second'},1,length(last450SplitCon)) repmat({'Second'},1,length(last450SplitExp))];
    groupCondition = [repmat({'Private'},1,length(first450SplitCon)) repmat({'Public'},1,length(first450SplitExp))...
         repmat({'Private'},1,length(last450SplitCon)) repmat({'Public'},1,length(last450SplitExp))];
    [p, tbl, stats, terms] = anovan(anovaSplit, {groupFirst, groupCondition}, 'model', 'interaction', 'varnames', {'First', 'Condition'});

    figure()
    subplot(1,2,1)
    a = first450SplitExp;
    b = last450SplitExp;
    [t, mean_a, mean_b, se_a, se_b] = ttest2_barplot(a,b,f);
    ylabel('Conflict')
    xlabel('Public')
    ylim([-100 100])
    subplot(1,2,2)
    a = first450SplitCon;
    b = last450SplitCon;
    [t, mean_a, mean_b, se_a, se_b] = ttest2_barplot(a,b,f);
    ylabel('Conflict')
    xlabel('Private')
    ylim([-100 100])
end

%% t-test between conditions for each timepoint
for ch=1:chunks
    [~,chunkPvals(ch,1),~,chunkStats(ch,:)] =ttest2(controlTC(ch,:),expTC(ch,:),'Vartype','unequal');  
end
clear ch

%% Multiple comparisons correction
% If we use Bonferroni here it is overly strict 0.05/chunks=sig. cutoff
% bfCorr=0.05/chunks;
% sigMask=chunkPvals(:,1)<=bfCorr; %Nothing is sig.

% If we use the classic FDR correction
% [h, crit_p, adj_ci_cvrg, adj_p]=fdr_bh(chunkPvals,cutoff,'pdep','yes');
[h crit_p]=fdr_bky(chunkPvals,cutoff,'yes');

%To find significant chunks of time
for r=1:chunks
    meanExpTP(r,1)= nanmean(expTC(r,:));
    meanConTP(r,1)=nanmean(controlTC(r,:));
end
clear r
fdrMask=find(h==1); %Timepoints (i.e., rows that are sig.)

%% Rearrange dataframe by group, rather than by participant ("halve" the data)
data1 = data(1:61,:);
data2 = data(62:122,:);

groupData = join(data1,data2, 'Keys', 'group');

groupData.Opp_conflict_total = (groupData.Opp_conflict_data1 + groupData.Opp_conflict_data2)/2;
groupData.conflict_total = (groupData.conflict_data1 + groupData.conflict_data2)/2;

% Consolidate conflict data as wanted
conflictData = table(groupData.group, groupData.experimental_data1,...
    groupData.Opp_conflict_data1, groupData.Opp_conflict_data2, groupData.conflict_data1, groupData.conflict_data2,...
    groupData.Opp_conflict_total, groupData.conflict_total, groupData.conflictBoth_data1,...
    groupData.carmaMeans_data1);
conflictData.Properties.VariableNames = {'group', 'experimental',...
    'selfReport1', 'selfReport2', 'rater1', 'rater2', 'selfReportAvg', 'raterAvg', 'raterOverall', 'carmaMeans'};

%% Sort by highest conflict convos
rankConflict = table(conflictData.group, conflictData.experimental);
rankConflict.Properties.VariableNames = {'group', 'experimental'};
rr = (1:height(rankConflict))';

% Sort by self-report avg
[conflictDataSorted_selfReport, rankSelfReportAvg] = sortrows(conflictData, [7,8], 'descend');
rr(rankSelfReportAvg) = rr;
rankConflict.rankSelfReportAvg = rr;
rr = (1:height(rankConflict))';
% Sort by rater avg
[conflictDataSorted_rater, rankRaterAvg] = sortrows(conflictData, [8,7], 'descend');
rr(rankRaterAvg) = rr;
rankConflict.rankRaterAvg = rr;
rr = (1:height(rankConflict))';
% Sort by rater overall
[conflictDataSorted_raterBoth, rankRaterBoth] = sortrows(conflictData, [9,10], 'descend');
rr(rankRaterBoth) = rr;
rankConflict.rankRaterBoth = rr;
rr = (1:height(rankConflict))';
% Sort by carma
[conflictDataSorted_carmaMeans, rankCarmaMeans] = sortrows(conflictData, [10,9], 'descend');
rr(rankCarmaMeans) = rr;
rankConflict.rankCarmaMeans = rr;
rr = (1:height(rankConflict))';

% Average over all 4 ranked methods to determine "average rank"
% Lower numbers are higher conflict convos
rankConflict.totalRank = (rankConflict.rankSelfReportAvg + rankConflict.rankRaterAvg + ...
    rankConflict.rankRaterBoth + rankConflict.rankCarmaMeans)/4;

%% Choose which average you want conflict sorting to be based on 
if conflictSort==1
    rankConflict = sortrows(rankConflict, [7]); % Over total average rank
elseif conflictSort==2  
    rankConflict = sortrows(rankConflict, [3]); % Over self report ranks
elseif conflictSort==3    
    rankConflict = sortrows(rankConflict, [4]); % Over rater average rank
elseif conflictSort==4
    rankConflict = sortrows(rankConflict, [6]); % Over CARMA average rank
elseif conflictSort==5
    rankConflict = sortrows(rankConflict, [5]); % Over rater average rank (overall)
end

%% Plot top conflict CARMA timecourses
% Note to self - Convo 40 is wonky
num = 5;
topConflictConvos = rankConflict.group(1:num);
topConflictExps = rankConflict.experimental(1:num);
botConflictConvos = rankConflict.group(62-num:61);
botConflictExps = rankConflict.experimental(62-num:61);

% topConflictConvos = {'58B','40'};
% topConflictExps = [0,1];
carmaIndexer = controlExp.group;

% figure()
% hold on
curTimecourses = nan(1163,num*2);
for i = 1:length(topConflictConvos)
    curGroupHigh{i} = topConflictConvos{i}; % group to be plotted
    curExpHigh{i} = topConflictExps(i); % exp or control
    if curExpHigh{i} == 1
        allExp{i} = 'High (Public)';
    else
        allExp{i} = 'High (Private)';
    end
    idx = find(strcmp(carmaIndexer,curGroupHigh{i})); % find appropriate carma timecourse index
    curTimecourse = carmaTimecourses{idx};
    curTimecourses(1:length(carmaTimecourses{idx}),i) = carmaTimecourses{idx};
%     if curExpHigh{i} == 1
%         p1 = plot(curTimecourse,'-','Color', [1 0 0]);
%     elseif curExpHigh{i}==0
%         p2 = plot(curTimecourse,'-','Color', [0 0 1]);
%     end
end
for i = 1:length(botConflictConvos)
    curGroupLow{i} = botConflictConvos{i}; % group to be plotted
    curExpLow{i} = botConflictExps(i); % exp or control
    if curExpLow{i} == 1
        allExp{i+num} = 'Low (Public)';
    else
        allExp{i+num} = 'Low (Private)';
    end
    idx = find(strcmp(carmaIndexer,curGroupLow{i})); % find appropriate carma timecourse index
    curTimecourse = carmaTimecourses{idx};
    curTimecourses(1:length(carmaTimecourses{idx}),i+num) = carmaTimecourses{idx};
%     if curExpLow{i} == 1
%         p3 = plot(curTimecourse,'-','Color', [1 .5 .5]);
%     elseif curExpLow{i}==0
%         p4 = plot(curTimecourse,'-','Color', [.5 .5 1]);
%     end
end

%% Plots for paper
% Plot Time course comparison of conditions
figure()
p1 = stdshade(expTC',.5,'r','r');
hold on
p2 = stdshade(controlTC',.5,'b','b');
ylim1 = -45;
ylim2 = 5;
ylim([ylim1,ylim2]);
annotation('doublearrow',[.075 .075],[.1,.93], 'LineWidth',1);
for i = 1:(length(fdrMask)-1) % shade in significant areas (4 x/y coordinates)
    fill([fdrMask(i) fdrMask(i)+1 fdrMask(i)+1 fdrMask(i)],[ylim2, ylim2, ylim1, ylim1], 'k','FaceAlpha', .3,'linestyle','none');
end
xline(450,'--k','LineWidth', 2)
legend([p1 p2],'Public','Private','Location','northwest')
ylabel('(Low)               Conflict               (High)','fontsize',14,'Position',[-70 -20 -1])
xlabel('Time (s)','fontsize',14)

% Plot Highest and lowest conflict convos
allGroups = [curGroupHigh curGroupLow];
cmap = hsv(num*2);
cmap2 = cmap([1:3,9:10, 4:8],:);

figure()
hold on
for i = 1:size(curTimecourses,2)
    pr{i} = plot(curTimecourses(1:900,i),'-', 'Color', cmap(i,:), 'LineWidth', .85);
%     pr{i} = plot(curTimecourses(1:900,i),'-', 'Color', cmap2(i,:), 'LineWidth', .8);
end
ylim([-100,100]);
ylabel('(Low)               Conflict               (High)','fontsize',14,'Position',[-88 -.000095 -1])
xlabel('Time (s)','fontsize',14)
yline(0,'--'); % Plot line at 0
annotation('doublearrow',[.06 .06],[.1,.93], 'LineWidth',1);
% legend([pr{1} pr{2} pr{3} pr{4} pr{5} pr{6} pr{7} pr{8} pr{9} pr{10}],...
%     allGroups, 'Orientation', 'horizontal', 'NumColumns',5)
legend([pr{1} pr{2} pr{3} pr{4} pr{5} pr{6} pr{7} pr{8} pr{9} pr{10}],...
    allExp, 'NumColumns',2,'fontsize',8,'Position',[0.49,0.75,0.41,0.18])

