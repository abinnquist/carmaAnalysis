%% Carma video analysis
% Will output statistics from the paper in the command window as well as
% create the figure 4 from the "Zoom Solution" paper.
% Make sure to run from the 'carmaAnalysis' folder that has helper scripts
% and data (carmaMeans & IPCdata_complete) stored in it. 

%% Load in external data 
carmaPath=cd; 
load(strcat(carmaPath,filesep,'carmaConvos.mat'))

trimvid=900; %Specified time at which to trim for analysis (some conversations ran longer).
cutoff=0.05; %FDR cutoff
%% Compile and average carma data
timeCourses=nan(1163,61); %Needs to be the length of the longest video or analysis will be incorrect
for v=1:61
    lenVid=length(carmaTimecourses{1,v});
    timeCourses(1:lenVid,v)=carmaTimecourses{1:lenVid,v};
end

%trim the timecourses  
timeChunk=timeCourses(1:trimvid,:); 
controlTC=timeChunk(:,controlExp.experimental(:)==0);
expTC=timeChunk(:,controlExp.experimental(:)==1);

%% t-test between conditions for each timepoint
for ch=1:trimvid
    [~,chunkPvals(ch,1),~,chunkStats(ch,:)] =ttest2(controlTC(ch,:),expTC(ch,:),'Vartype','unequal');  
end
clear ch

%% Multiple comparisons correction
% you can change the output to get the critical p-val if wanted
% FDR correction (crit_p=0.0182)
[h, ~]=fdr_bky(chunkPvals,cutoff,'yes');

%To find significant chunks of time
for r=1:trimvid
    meanExpTP(r,1)= nanmean(expTC(r,:));
    meanConTP(r,1)=nanmean(controlTC(r,:));
end
clear r
fdrMask=find(h==1); %Timepoints (i.e., rows that are sig.)

%Outputs the first timepoint for a sig. period of time longer than 10 seconds
count=1;
for i=1:length(fdrMask)
    if diff(fdrMask(i:i+1))==1
        count=count+1;
        if count > 10
            startCh=fdrMask(i-9,1);
            pval=chunkPvals(startCh);
            fprintf('Significant chunk of time starting at time: %d(s), p = %.3f \n',startCh,pval)
            break
        end
    else
        count=1;
    end
end
clearvars count v i lenVid cutoff

%% t-Test between condition X first/second half  
% Will output the t-stat & p-val for the first and second half of the
% conversations for difference of condition.
a = nanmean(controlTC(1:450,:)); b = nanmean(expTC(1:450,:));
c = nanmean(controlTC(451:900,:)); d = nanmean(expTC(451:900,:));

[~,pvalsHalf(1,1),~,halfStats(1,:)] =ttest2(a,b,'Vartype','unequal'); 
[~,pvalsHalf(2,1),~,halfStats(2,:)] =ttest2(c,d,'Vartype','unequal'); 

%To get Cohen's D
for t=1:2
    n1=30; n2=31;
    if t==1
        m1=nanmean(a); m2=nanmean(b);s1=nanvar(a); s2=nanvar(b);   
    else
        m1=nanmean(c); m2=nanmean(d); s1=nanvar(c); s2=nanvar(d);   
    end
    pooledSD = sqrt((((n1-1)*s1) + ((n2-1)*s2))/(n1 + n2 - 2));
    cohendHalf(t) = (m1-m2) / pooledSD;  
end

fprintf('Conflict in the public condition vs the private condition, first half (t = %.3f, p = %.3f, d = %.3f)\n',...
    halfStats(1).tstat,pvalsHalf(1),cohendHalf(1))

fprintf('Conflict in the public condition vs the private condition, second half (t = %.3f, p = %.3f, d = %.3f)\n',...
    halfStats(2).tstat,pvalsHalf(2),cohendHalf(2))

%% Average CARMA rating difference between conditions
carmaMeans=nanmean(timeCourses(1:1163,:))';
a = carmaMeans(controlExp.experimental == 1);
b = carmaMeans(controlExp.experimental == 0);
[~, p, ~, stats] = ttest2(a, b);

%To get Cohen's D
n1=30; n2=31;
m1=nanmean(a); m2=nanmean(b);s1=nanvar(a); s2=nanvar(b);   
pooledSD = sqrt((((n1-1)*s1) + ((n2-1)*s2))/(n1 + n2 - 2));
cohend = (m1-m2) / pooledSD;  

fprintf('More conflict in the public condition vs the private condition (t = %.3f, p = %.3f, d = %.3f)\n',...
    stats.tstat,p,cohend)

%% Rearrange dataframe by group, rather than by participant ("halve" the data)
opts = detectImportOptions('IPCdata_complete.csv');
opts = setvartype(opts, 'group', 'char');  %or 'string' if you prefer
data = readtable('IPCdata_complete.csv', opts);

isExp = logical(data.experimental);
expData = data(isExp,:);
controlData = data(~isExp,:);

data1 = data(1:61,:);
data2 = data(62:122,:);

groupData = join(data1,data2, 'Keys', 'group');

% Consolidate conflict data as wanted
conflictData = table(groupData.group, groupData.experimental_data1, groupData.conflictBoth_data1,...
    groupData.carmaMeans_data1);
conflictData.Properties.VariableNames = {'group', 'experimental','raterOverall', 'carmaMeans'};

%% Sort by highest conflict convos
rankConflict = table(conflictData.group, conflictData.experimental);
rankConflict.Properties.VariableNames = {'group', 'experimental'};
rr = (1:height(rankConflict))';

% Sort by carma
[conflictDataSorted_carmaMeans, rankCarmaMeans] = sortrows(conflictData, [4,3], 'descend');
rr(rankCarmaMeans) = rr;
rankConflict.rankCarmaMeans = rr;
rr = (1:height(rankConflict))';
rankConflict = sortrows(rankConflict, 3); % Over CARMA average rank

%% Plot top conflict CARMA timecourses
num = 5;
topConflictConvos = rankConflict.group(1:num);
topConflictExps = rankConflict.experimental(1:num);
botConflictConvos = rankConflict.group(62-num:61);
botConflictExps = rankConflict.experimental(62-num:61);

carmaIndexer = controlExp.group;

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
ylabel('Friendly               Neutral               Conflict','fontsize',14,'Position',[-70 -20 -1])
xlabel('Time (s)','fontsize',14)

% Plot Highest and lowest conflict convos
allGroups = [curGroupHigh curGroupLow];
cmap = hsv(num*2);
cmap2 = cmap([1:3,9:10, 4:8],:);

figure()
hold on
for i = 1:size(curTimecourses,2)
    pr{i} = plot(curTimecourses(1:900,i),'-', 'Color', cmap(i,:), 'LineWidth', .85);
end
ylim([-100,100]);
ylabel('Friendly               Neutral               Conflict','fontsize',14,'Position',[-88 -.000095 -1])
xlabel('Time (s)','fontsize',14)
yline(0,'--'); % Plot line at 0
annotation('doublearrow',[.06 .06],[.1,.93], 'LineWidth',1);
legend([pr{1} pr{2} pr{3} pr{4} pr{5} pr{6} pr{7} pr{8} pr{9} pr{10}],...
    allExp, 'NumColumns',2,'fontsize',8,'Position',[0.49,0.75,0.41,0.18])

