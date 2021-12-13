%Plot figures for Scenario 1 - Case 9
%Author: Sebastian A. Nugroho
%Date: 3/14/2021

clear
close all

%Load data from simulation
load('Data_scenario1_Case9_dist_low.mat');

%Assign values
t1_l = t2; %no control
t2_l = t4; %NDAE
t3_l = t5; %LQR
t4_l = t6; %AGC
x1_l = x_pre2; %no control
x2_l = x_res2; %NDAE
x3_l = x_res3; %LQR
x4_l = x_res4; %AGC

%Load data from simulation
load('Data_scenario1_Case9_dist_moderate.mat');

%Assign values
t1_m = t2; %no control
t2_m = t4; %NDAE
t3_m = t5; %LQR
t4_m = t6; %AGC
x1_m = x_pre2; %no control
x2_m = x_res2; %NDAE
x3_m = x_res3; %LQR
x4_m = x_res4; %AGC

%Load data from simulation
load('Data_scenario1_Case9_dist_high.mat');

%Assign values
t1_h = t2; %no control
t2_h = t4; %NDAE
t3_h = t5; %LQR
t4_h = t6; %AGC
x1_h = x_pre2; %no control
x2_h = x_res2; %NDAE
x3_h = x_res3; %LQR
x4_h = x_res4; %AGC

%Compute frequency deviation
stateidx = 4; %Generator 1 frequency
t_10 = 10;
idx_10 = find(t4 == t_10);
t_20 = 20;
idx_20 = find(t4 == t_20);
t_30 = 30;
idx_30 = find(t4 == t_30);

err_low_NDAE = norm(case_sys.w0*ones(1,case_sys.nG)-x2_l(idx_10,stateidx:stateidx+2),2);
err_low_LQR = norm(case_sys.w0*ones(1,case_sys.nG)-x3_l(idx_10,stateidx:stateidx+2),2);
err_low_AGC = norm(case_sys.w0*ones(1,case_sys.nG)-x4_l(idx_10,stateidx:stateidx+2),2);
err_low_NDAE = [err_low_NDAE norm(case_sys.w0*ones(1,case_sys.nG)-x2_l(idx_20,stateidx:stateidx+2),2)];
err_low_LQR = [err_low_LQR norm(case_sys.w0*ones(1,case_sys.nG)-x3_l(idx_20,stateidx:stateidx+2),2)];
err_low_AGC = [err_low_AGC norm(case_sys.w0*ones(1,case_sys.nG)-x4_l(idx_20,stateidx:stateidx+2),2)];
err_low_NDAE = [err_low_NDAE norm(case_sys.w0*ones(1,case_sys.nG)-x2_l(idx_30,stateidx:stateidx+2),2)];
err_low_LQR = [err_low_LQR norm(case_sys.w0*ones(1,case_sys.nG)-x3_l(idx_30,stateidx:stateidx+2),2)];
err_low_AGC = [err_low_AGC norm(case_sys.w0*ones(1,case_sys.nG)-x4_l(idx_30,stateidx:stateidx+2),2)];

err_mod_NDAE = norm(case_sys.w0*ones(1,case_sys.nG)-x2_m(idx_10,stateidx:stateidx+2),2);
err_mod_LQR = norm(case_sys.w0*ones(1,case_sys.nG)-x3_m(idx_10,stateidx:stateidx+2),2);
err_mod_AGC = norm(case_sys.w0*ones(1,case_sys.nG)-x4_m(idx_10,stateidx:stateidx+2),2);
err_mod_NDAE = [err_mod_NDAE norm(case_sys.w0*ones(1,case_sys.nG)-x2_m(idx_20,stateidx:stateidx+2),2)];
err_mod_LQR = [err_mod_LQR norm(case_sys.w0*ones(1,case_sys.nG)-x3_m(idx_20,stateidx:stateidx+2),2)];
%err_mod_AGC = [err_mod_AGC norm(case_sys.w0*ones(1,case_sys.nG)-x4_m(idx_20,stateidx:stateidx+2),2)];
err_mod_NDAE = [err_mod_NDAE norm(case_sys.w0*ones(1,case_sys.nG)-x2_m(idx_30,stateidx:stateidx+2),2)];
err_mod_LQR = [err_mod_LQR norm(case_sys.w0*ones(1,case_sys.nG)-x3_m(idx_30,stateidx:stateidx+2),2)];
%err_mod_AGC = [err_mod_AGC norm(case_sys.w0*ones(1,case_sys.nG)-x4_m(idx_30,stateidx:stateidx+2),2)];

err_high_NDAE = norm(case_sys.w0*ones(1,case_sys.nG)-x2_h(idx_10,stateidx:stateidx+2),2);
%err_high_LQR = norm(case_sys.w0*ones(1,case_sys.nG)-x3_h(idx_10,stateidx:stateidx+2),2);
%err_mod_AGC = norm(case_sys.w0*ones(1,case_sys.nG)-x4_h(idx_10,stateidx:stateidx+2),2);
err_high_NDAE = [err_high_NDAE norm(case_sys.w0*ones(1,case_sys.nG)-x2_h(idx_20,stateidx:stateidx+2),2)];
%err_high_LQR = [err_high_LQR norm(case_sys.w0*ones(1,case_sys.nG)-x3_h(idx_20,stateidx:stateidx+2),2)];
%err_mod_AGC = [err_mod_AGC norm(case_sys.w0*ones(1,case_sys.nG)-x4_h(idx_20,stateidx:stateidx+2),2)];
err_high_NDAE = [err_high_NDAE norm(case_sys.w0*ones(1,case_sys.nG)-x2_h(idx_30,stateidx:stateidx+2),2)];
%err_high_LQR = [err_high_LQR norm(case_sys.w0*ones(1,case_sys.nG)-x3_h(idx_30,stateidx:stateidx+2),2)];
%err_mod_AGC = [err_mod_AGC norm(case_sys.w0*ones(1,case_sys.nG)-x4_h(idx_30,stateidx:stateidx+2),2)];


%Plot figures
stateidx = 4; %Generator 1 frequency
h(1) = figure;
fs = 14;
%Plot time index
t_final = 55;
idx_fin = find(t4 == t_final);
hold on
p(1) = plot(t2_l(1:idx_fin),x2_l(1:idx_fin,stateidx)/(2*pi),'LineWidth',1.2,'Color',[0.88 0.00 0.15]);
p(2) = plot(t3_l(1:idx_fin),x3_l(1:idx_fin,stateidx)/(2*pi),'LineWidth',1.2,'Color',[0.47 0.67 0.19],'LineStyle','-.');
p(3) = plot(t4_l(1:idx_fin),x4_l(1:idx_fin,stateidx)/(2*pi),'LineWidth',1.2,'Color',[0.00 0.45 0.74],'LineStyle','--');
set(gcf,'color','w');
title('\rm 9-Bus Network','FontSize',fs-5)
ax = gca;
ax.TitleFontSizeMultiplier = 1;
hold off
grid on
box on
xlim([0 t2_l(idx_fin)]);
ylim([59.9975 60.0001]);
xlabel('$t$ [sec]', 'interpreter','latex','FontName','Times New Roman','FontSize',fs);
ylabel('$\frac{1}{2\pi}\omega$ [Hz]', 'interpreter','latex','FontName','Times New Roman','FontSize',fs);
set(gca,'fontsize',fs-3);
leg1 = legend(p,{'NDAE-control','LQR-control','AGC'});
set(leg1,'Interpreter','latex');
set(leg1,'FontSize',fs-4);
set(leg1,'Location','southeast');
%set(leg1,'Orientation','horizontal');

%Add text
str = '$\rho_{\mathrm{L}} = 0.04$';
text(0.022727272727273,0.940238344579429,str,'Interpreter','latex','FontSize',fs-3,'Units','normalized');

set(h(1), 'Position', [300 100 390 240])
print(h(1), 'scenario1_case9_low.eps', '-depsc2','-r600')

h(2) = figure;
%Plot time index
t_final = 55;
idx_fin = find(t4 == t_final);
hold on
p(1) = plot(t2_m(1:idx_fin),x2_m(1:idx_fin,stateidx)/(2*pi),'LineWidth',1.2,'Color',[0.88 0.00 0.15]);
p(2) = plot(t3_m(1:idx_fin),x3_m(1:idx_fin,stateidx)/(2*pi),'LineWidth',1.2,'Color',[0.47 0.67 0.19],'LineStyle','--');
p(3) = plot(t4_m,x4_m(:,stateidx)/(2*pi),'LineWidth',1.2,'Color',[0.00 0.45 0.74],'LineStyle','-.');
set(gcf,'color','w');
title('\rm 9-Bus Network','FontSize',fs-5)
ax = gca;
ax.TitleFontSizeMultiplier = 1;
hold off
grid on
box on
xlim([0 t2_m(idx_fin)]);
ylim([59.994 60.0001]);
xlabel('$t$ [sec]', 'interpreter','latex','FontName','Times New Roman','FontSize',fs);
ylabel('$\frac{1}{2\pi}\omega$ [Hz]', 'interpreter','latex','FontName','Times New Roman','FontSize',fs);
set(gca,'fontsize',fs-3);
leg1 = legend(p,{'NDAE-control','LQR-control','AGC'});
set(leg1,'Interpreter','latex');
set(leg1,'FontSize',fs-4);
set(leg1,'Location','southeast');
%set(leg1,'Orientation','horizontal');

%Add text
str = '$\rho_{\mathrm{L}} = 0.08$';
text(0.022727272727273,0.940238344579429,str,'Interpreter','latex','FontSize',fs-3,'Units','normalized');

set(h(2), 'Position', [300 100 390 240])
print(h(2), 'scenario1_case9_mod.eps', '-depsc2','-r600')

h(3) = figure;
%Plot time index
t_final = 55;
idx_fin = find(t4 == t_final);
hold on
p(1) = plot(t2_h(1:idx_fin),x2_h(1:idx_fin,stateidx)/(2*pi),'LineWidth',1.2,'Color',[0.88 0.00 0.15]);
p(2) = plot(t3_h,x3_h(:,stateidx)/(2*pi),'LineWidth',1.2,'Color',[0.47 0.67 0.19],'LineStyle','-.');
p(3) = plot(t4_h,x4_h(:,stateidx)/(2*pi),'LineWidth',1.2,'Color',[0.00 0.45 0.74],'LineStyle','--');
set(gcf,'color','w');
title('\rm 9-Bus Network','FontSize',fs-5)
ax = gca;
ax.TitleFontSizeMultiplier = 1;
hold off
grid on
box on
xlim([0 t2_h(idx_fin)]);
ylim([59.993 60.0001]);
xlabel('$t$ [sec]', 'interpreter','latex','FontName','Times New Roman','FontSize',fs);
ylabel('$\frac{1}{2\pi}\omega$ [Hz]', 'interpreter','latex','FontName','Times New Roman','FontSize',fs);
set(gca,'fontsize',fs-3);
leg1 = legend(p,{'NDAE-control','LQR-control','AGC'});
set(leg1,'Interpreter','latex');
set(leg1,'FontSize',fs-4);
set(leg1,'Location','southeast');
%set(leg1,'Orientation','horizontal');

%Add text
str = '$\rho_{\mathrm{L}} = 0.12$';
text(0.022727272727273,0.940238344579429,str,'Interpreter','latex','FontSize',fs-3,'Units','normalized');

set(h(3), 'Position', [300 100 390 240])
print(h(3), 'scenario1_case9_high.eps', '-depsc2','-r600')