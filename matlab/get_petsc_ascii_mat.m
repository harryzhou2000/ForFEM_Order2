

fin = fopen('mat1.txt');

siz = 0;
is = nan(100,1);
js = nan(100,1);
vs = nan(100,1);
while(~feof(fin))
    line = fgetl(fin);
    row = regexp(line,'row *(\d+):','tokens');
    row = str2double(row{1})+1;
    pairs = regexp(line,'\((\d+), *([+-]?\d*\.?\d+[eE]?[+-]?\d*)\)','tokens');
    for i = 1:numel(pairs)
        col = str2double(pairs{i}{1})+1;
        val = str2double(pairs{i}{2});
        siz = siz + 1;
        js(siz) = col;
        is(siz) = row;
        vs(siz) = val;
    end
end


fclose(fin);

mat = sparse(is,js,vs);
spy(mat)