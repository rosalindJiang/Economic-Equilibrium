% generate pi and c

clc, clear all, close all

m = 100;
n = 1000;

% the members in pi follows the standard normal distribution
pi = randn([n,m+1]);

% the members in c follows the gamma distribution
% change the variance of c
c_1over2 = gamrnd(2,(1/2),[1,n]);  % mean=1, var=1/2
c_1over3 = gamrnd(3,(1/3),[1,n]);  % mean=1, var=1/3
c_1over5 = gamrnd(5,(1/5),[1,n]);  % mean=1, var=1/5
c_1over9 = gamrnd(9,(1/9),[1,n]);  % mean=1, var=1/9

% save the data
save("generator.mat",'pi','c_1over2','c_1over3','c_1over5','c_1over9');