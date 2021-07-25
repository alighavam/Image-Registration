function [patFlag,pat] = isPat(pat)

pat = str2num(pat(4:end));
if ((pat>=11 && pat<=21))
    patFlag = 1;
else
    patFlag = 0;
end