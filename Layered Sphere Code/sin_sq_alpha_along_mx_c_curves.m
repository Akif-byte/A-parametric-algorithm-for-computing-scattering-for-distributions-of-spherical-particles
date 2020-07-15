clc
clear all

%This script calculates the Average periods, phase of alpha/beta for a given c
%It calls the function calc_period which takes in the inputs c and the mode
%number n
%Range of c which will be c_min=min(m1)*min(x), c_max=max(m1)*max(x)
%Range of m1 is chosen to be 1.2<m1<1.4
c=[70];
m1=[1.2:0.0001:1.4];
%n is the mode number
n=1;
for i=1:length(c)
    [Avg_period_alpha_mode_1(i),Avg_period_beta_mode_1(i),Phase_alpha_mode_1(i),Phase_beta_mode_1(i)]=calc_period(c(i),m1,n);
end


