function auto_correlation_plot(x)

N = length(x);
res_cov = xcorr(x-mean(x),'coeff');
[max_res_cov, start_i]=max(res_cov);
s_cov = max_res_cov/sqrt(N);
% s_cov = std(res_cov);

figure('position',[0,0,350,250])
plot(res_cov(start_i:end),':','linewidth',2); hold on; grid on;
plot(2*s_cov*ones(N,1),'r--','linewidth',2);
plot(-2*s_cov*ones(N,1),'r--','linewidth',2);
ylabel('R{vv}'); xlabel('lag index');
legend('autocorrelation estimates','2-sigma confidence interval');
% ylim([-30*s_cov,80*s_cov]);
end