function [t, mean_a, mean_b, se_a, se_b] = ttest2_barplot(a,b,c)
    se_a = stderr(a);
    se_b = stderr(b);
    mean_a = mean(a,'omitnan');
    mean_b = mean(b,'omitnan');
    c1 = c(1);
    c2 = c(2);
    
    means = [mean_a, mean_b];
    err = [se_a, se_b];
    
    b = bar(c1,mean_a);
    b(1).FaceColor = [.5,1,.7];
%     b(1).FaceColor = [.5,.7,1];

    hold on
    b2 = bar(c2,mean_b);
    b2(1).FaceColor = [1,.5,.7];
%     b2(1).FaceColor = [1,.7,.5];



    
    t = errorbar(c,means, err, err);
    t.Color = [0 0 0];                            
    t.LineStyle = 'none';  
    hold off
end