
function TCtable = writecarmacsv(timecourse,conExp)
    tNames=[];
    for i=1:length(timecourse)
        num=num2str(i);
        vidNum=strcat("t_",num);
        tNames=[tNames;vidNum];
    end

    tNames=tNames';

    TCtable=array2table(timecourse');
    TCtable.Properties.VariableNames=tNames;
    TCtable=[TCtable,conExp];
end