function [cohensd] = computepairedcohensd(premetric,postmetric)

paireddiffs = postmetric - premetric;

simpled = nanmean(paireddiffs)/nanstd(paireddiffs);

r = corrcoef(premetric, postmetric,'rows','pairwise');
r = r(2);

cohensd = simpled/sqrt(1-r);