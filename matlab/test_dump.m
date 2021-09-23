% A=f_get_petsc_ascii_mat('cooler_elas_S_A.dump');
% 
% M=f_get_petsc_ascii_mat('cooler_elas_S_M.dump');
    
% A = PetscBinaryRead('cooler_elas_S_A.dump');
% M = PetscBinaryRead('cooler_elas_S_M.dump');
A = PetscBinaryRead('cooler_elas_origin_A.dump');
M = PetscBinaryRead('cooler_elas_origin_M.dump');
%%
M(M>0&M<1e-50)=1;
%%
dA = diag(A);
[vminda,iminda] = min(dA)
dM = diag(M);
% dM(dM>1e50) = 1e-50;
%%
tic;    
[v,d] = eigs(A,M,10,'smallestabs','Tolerance',1e-1);
fprintf('solve time: %g\n',toc);
%%
sM = sum(M,1);
sM(sM>1e50) = 1e-50;
M2 = diag(sM); max(sM)
tic;
[v,d] = eigs(A,M2,10,'smallestabs','Tolerance',1e-6);
fprintf('solve time: %g\n',toc);

%%
[v_m,d_m] = eigs(M,10,'smallestabs');
%%
[v_k,d_k] = eigs(A,10,'smallestabs');
%%
v1 = reshape(v(:,1),3,[]);
U1 = v1(1,:);
V1 = v1(2,:);
W1 = v1(3,:);
%%
fout = fopen('modeExport.txt','w');
fprintf(fout,'%.16e,%.16e,%.16e\n',v1);
fclose(fout);
