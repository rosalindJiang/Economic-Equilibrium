% validation

clc,clear all,close all

% Gdt1 = cell2mat(struct2cell(load("PE1_c_var_1over3.mat","G")));
% Gdt2 = cell2mat(struct2cell(load("PE1_c_var_1over3dt2.mat","G")));
% 
% time1 = 1:1:100000;
% time2 = 2:2:100000;
% 
% x = time1;
% y1 = Gdt1;
% p1 = plot(x,y1,"--",'DisplayName','dt=0.0001');
% p1.LineWidth = 1.5; 
% hold on;
% 
% x = time2;
% y2 = Gdt2;
% p2 = plot(x,y2,":",'DisplayName','dt=0.0002');
% p2.LineWidth = 1.5; 
% hold on;
% 
% 
% legend('FontSize',17)
% 
% title('The Economic Growth Rate with Variance of c 1/3 for PE1 for Different Small Time Step')
% xlabel('Time') 
% ylabel('Economic Growth Rate') 
% ylim([0 4])
% ax = gca;
% ax.FontSize = 20;


ej_pe1 = cell2mat(struct2cell(load("PE1_c_var_1over3.mat","e")));
ej_pe2 = cell2mat(struct2cell(load("PE2_c_var_1over3.mat","e")));

time = 1:1:100000;

x = time;

subplot(2,1,1)
y1 = ej_pe1';
p1 = plot(x,y1);
% p1.LineWidth = 1.5; 
xlabel('Time') 
ylabel('Excess Demand for Non-Capital Good') 
title('The Excess Demand for Non-Capital Good with Variance of c 1/3 under PE1')
ax = gca;
ax.FontSize = 20;
hold on;

subplot(2,1,2)
y2 = ej_pe2';
p2 = plot(x,y2);
% p2.LineWidth = 1.5; 
xlabel('Time') 
ylabel('Excess Demand for Non-Capital Good') 
title('The Excess Demand for Non-Capital Good with Variance of c 1/3 under PE2')
ylim([-0.1,0.1])
ax = gca;
ax.FontSize = 20;
hold on;



