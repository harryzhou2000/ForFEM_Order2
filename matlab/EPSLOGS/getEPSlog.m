function log = getEPSlog(filename)

pdouble = '[+-]?\d*\.?\d+[eE]?[+-]?\d*';
pint = '[+-]?\d+';
log = nan(1,3);
nlog = 0;

fin = fopen(filename);
while(~feof(fin))
    line = fgetl(fin);

    nums = regexp(line,['(',pdouble,')'],'tokens');
    if(numel(nums)==4)
        nlog = nlog+1;
        log(nlog,1) = str2double(nums{2});
        log(nlog,2) = str2double(nums{3});
        log(nlog,3) = str2double(nums{4});
    end
end
