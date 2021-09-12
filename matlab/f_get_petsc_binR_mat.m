function mat = f_get_petsc_binR_mat(path)

fin = fopen(path,'r','b');
nrow = fread(fin, 1, 'int32');
rowsiz = fread(fin, rowsiz, 'int32');
i_csr = [0,cumsum(rowsiz)];
Petsc




fclose(fin);