%This function caluclates the coefficients pi &qi of the least squares fit of
%the the period T as a function of c for homogeneous spheres. The
%coefficients pi and qi are stored in the vectors Coeff_vec1 and
%Coeff_vec2 respectively. 
clc
clear all

%Range of c, m_min*x_min<c<m_max*x_max. For the sample distribution we have
%taken 1.2<m<1.8 and 10<x<20 
c=[10:0.001:40];
m=[1:0.001:1.8];
%n is the mode number
n=1;
for i=1:length(c)
    %calc_period calculates the Average periods of sin^2(alpha) and
    %sin^2(beta) for a given c and mode n. m is assumed to be in the
    %range of 1.2<m<1.8
    [Avg_period_alpha(i),Avg_period_beta(i)]=calc_period(c(i),n,m);
    Avg_period_alpha=Avg_period_alpha;
    Avg_period_beta=Avg_period_beta;
end
%This calculates the coefficients p1,p2,p3,p4 & q1,q2,q3,q4 and stores
%them in the vectors Coeff_vec1 & Coeff_vec2 respectively and saves
%them. So the pi,qi for the first mode is stored in the file
%Coeff_vec1_n_1, Coeff_vec2_n_1, for the second mode the file names are
%Coeff_vec1_n_2, Coeff_vec2_n_2 and so on for higher modes...
F1=[1./(2*c);1./(2*c.^2);sin(2*c+(n-1)*pi)./(2*c);sin(2*c+(n-1)*pi)./(2*c.^2)]';
F2=[1./(2*c);1./(2*c.^2);sin(2*c+n*pi)./(2*c);sin(2*c+n*pi)./(2*c.^2)]';
Coeff_vec1=F1\(Avg_period_alpha-pi)';
Coeff_vec2=F2\(Avg_period_alpha-pi)';
save(sprintf('Coeff_vec1_n_%d', n),'Coeff_vec1');
save(sprintf('Coeff_vec2_n_%d', n),'Coeff_vec2');






