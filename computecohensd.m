function [cohensd, cilow, cihigh] = computecohensd(nn1, mean1, std1, nn2, mean2, std2)
%n1 and n2 are group sizes here
n1 = nn1-1;
n2 = nn2-1;

stdpooled = sqrt((std1^2*n1+std2^2*n2)/(n1+n2));
cohensd  = abs((mean1-mean2)/stdpooled);



temp1 = (nn1+nn2)./(nn1*nn2);
temp2 = cohensd^2./(2*(n1+n2));
temp3 = (nn1+nn2)./(n1+n2);
stderror = (temp1 + temp2).*temp3

cilow = cohensd-stderror*1.96
cihigh = cohensd+stderror*1.96