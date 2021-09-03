
log_chol_elas=getKSPlog('BEAM_chol_elas_KSP.log');
log_icc_elas=getKSPlog('BEAM_icc_elas_KSP.log');

figure(1);clf;
hold on;
plot(log_chol_elas/log_chol_elas(1))
plot(log_icc_elas/log_chol_elas(1))
set(gca,'YScale','log');
xlabel('iteration');
title('unconverged relative residual');
legend('elastic cholesky','elastic icc');
grid on;
% 1.39813140906335 icc:1.31824821686314 

%%
log_chol_ther=getKSPlog('BEAM_chol_ther_KSP.log');
log_icc_ther=getKSPlog('BEAM_icc_ther_KSP.log');

figure(1);clf;
hold on;
plot(log_chol_ther/log_chol_ther(1))
plot(log_icc_ther/log_icc_ther(1))
set(gca,'YScale','log');
xlabel('iteration');
title('unconverged relative residual');
legend('thermal cholesky','thermal icc');
grid on;
% 15.6204882242091 icc:54.3898142039980

