function log = getKSPlog(filename)

pdouble = '[+-]?\d*\.?\d+[eE]?[+-]?\d*';
pint = '[+-]?\d+';
log = nan(1,1);
nlog = 0;

fin = fopen(filename);
while(~feof(fin))
    line = fgetl(fin);
    nums = regexp(line,['(',pdouble,')'],'tokens');
    if(numel(nums)==2)
        nlog = nlog+1;
        log(nlog,1) = str2double(nums{2});
    end
end
