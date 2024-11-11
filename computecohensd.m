function [cohensd] = computecohensd(n1, mean1, std1, n2, mean2, std2)
%n1 and n2 are group sizes here
n1 = n1-1;
n2 = n2-1;

stdpooled = sqrt((std1^2*n1+std2^2*n2)/(n1+n2));
cohensd  = (mean1-mean2)/stdpooled;