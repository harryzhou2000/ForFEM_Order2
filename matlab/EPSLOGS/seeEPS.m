
log_rqcg_elas=getEPSlog('BEAM_rqcg_elas_EPS.log');
log_lobpcg_elas=getEPSlog('BEAM_lobpcg_elas_EPS.log');

figure(1);clf;
subplot(2,1,1);cla;hold on;
plot(log_rqcg_elas(:,3))
plot(log_lobpcg_elas(:,3))
set(gca,'YScale','log');
xlabel('iteration');
title('unconverged relative residual');
legend('elastic rqcg','elastic lobpcg');
grid on;
subplot(2,1,2);cla;hold on;
plot(log_rqcg_elas(:,1))
plot(log_lobpcg_elas(:,1))
xlabel('iteration');
title('converged number of dimensions');
legend('elastic rqcg','elastic lobpcg');
grid on;
%%


log_rqcg_ther=getEPSlog('BEAM_rqcg_ther_EPS.log');
log_lobpcg_ther=getEPSlog('BEAM_lobpcg_ther_EPS.log');

figure(1);clf;
subplot(2,1,1);cla;hold on;
plot(log_rqcg_ther(:,3))
plot(log_lobpcg_ther(:,3))
set(gca,'YScale','log');
xlabel('iteration');
title('unconverged relative residual');
legend('thermal rqcg','thermal lobpcg');
grid on;
subplot(2,1,2);cla;hold on;
plot(log_rqcg_ther(:,1))
plot(log_lobpcg_ther(:,1))
xlabel('iteration');
title('converged number of dimensions');
legend('thermal rqcg','thermal lobpcg');
grid on;

