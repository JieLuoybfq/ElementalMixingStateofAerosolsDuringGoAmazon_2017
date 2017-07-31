function SimpleStats = SimpleStats(RowVector)
% SimpleStats = SimpleStats(RowVector)
% Calculates simple statistics on a row vector and outputs a
% structure of those statistics

t = 1.96; %used to find significance using students t-distribution @95% confidence
n = length(RowVector);
Max = max(RowVector);
Min = min(RowVector);
Avg = nanmean(RowVector);
Median = nanmedian(RowVector);
Mode = mode(RowVector);
Stddev = nanstd(RowVector);
Variance = nanvar(RowVector);
ConInt = t.*Stddev./sqrt(n);

SimpleStats = struct('Max',Max,'Min',Min,'Avg',Avg,'Median',Median,'Mode',Mode,'Stddev',Stddev,'Variance',Variance,'ConInt',ConInt,'number',n);
end

