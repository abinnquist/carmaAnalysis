function s = stderr(a) 
    t = size(a);
    if t(1) == 1
        l = t(2);
        s = std(a, 'omitnan')/sqrt(l);
    elseif t(2) == 1
        l = t(1);
        s = std(a, 'omitnan')/sqrt(l);
    else
        for i = 1:t(2)
            l(i) = sum(~isnan(a(:,i)));
            if l(i) == 0
                s(i) = Inf;
            else
                s(i) = std(a(:,i),'omitnan')/sqrt(l(i));
            end
        end
    end

end