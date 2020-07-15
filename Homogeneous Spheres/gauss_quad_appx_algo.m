%Gauss quadrature for evauluation of the integral I using the trigonometric approximation 
%a,b limits of outer integral
%c,d limts of inner integral
%p(m,x)=pdf_ij the probability density function
function [Int_eva]=gauss_quad_appx_algo(a,b,c,d,Avg_period_alpha_mode_1,Avg_period_beta_mode_1,Avg_period_alpha_mode_2,Avg_period_beta_mode_2,Avg_period_alpha_mode_3,Avg_period_beta_mode_3,wi,xi)
%wi,xi are Weights and nodes for Gauss quadrature we have used the same
%number of weights and nodes for both parameters m and x

% Converting limits of integration to -1 to 1
p1=(b-a)/2;
p2=(b+a)/2;
q1=(d-c)/2;
q2=(d+c)/2;
x=p1*xi+p2;
m=q1*xi+q2;
Int_eva=0;

%Uniform Distribution
pdf_ij=1;
%Normal Distribution Parameters
% mu_x=15; mu_m=1.5; sig_x=1.67; sig_m=0.1;

%Bi-Modal Distribution Parameters
% mu1_x=13;mu2_x=17; mu1_m=1.4;mu2_m=1.6; sig1_x=1;sig2_x=1; sig1_m=0.06;sig2_m=0.06; 
% p=0.5; %Weight for B-Modal distribution; f= p*Normal1 + (1-p)*Normal2

%Calculating values of Average period for a given c=m*x,
%Average_period_alpha/beta(1) corresponds to c=12,
%Average_period_alpha/beta(end) corresponds to c=36

c=m*x';
Index=floor((c-12)*1000)+1;
T1_mode_1=Avg_period_alpha_mode_1(Index);
T2_mode_1=Avg_period_beta_mode_1(Index);
T1_mode_2=Avg_period_alpha_mode_2(Index);
T2_mode_2=Avg_period_beta_mode_2(Index);
T1_mode_3=Avg_period_alpha_mode_3(Index);
T2_mode_3=Avg_period_beta_mode_3(Index);

%gaussian quadrature over the entire parametric space in m and x but using
%the approximate trigonometric function forms instead of the bessel
%function form of |a_n|^2 & |b_n|^2

for i=1:length(m)
    for j=1:length(x)
        
%       Normal distribution
%       pdf_ij=1/(2*pi*sig_x*sig_m)*exp(-0.5*((((x(j)-mu_x)/sig_x)^2) + ((m(i)-mu_m)/sig_m)^2));
        
%       Bi-Modal distribution
%         pdf1_ij=1/(2*pi*sig1_x*sig1_m)*exp(-0.5*((((x(j)-mu1_x)/sig1_x)^2) + ((m(i)-mu1_m)/sig1_m)^2)); 
%         pdf2_ij=1/(2*pi*sig2_x*sig2_m)*exp(-0.5*((((x(j)-mu2_x)/sig2_x)^2) + ((m(i)-mu2_m)/sig2_m)^2));
%         pdf_ij=p*pdf1_ij+(1-p)*pdf2_ij;
        
%       Multiplying by the weights of Gauss-Legendre quadrature wi
        Int_eva=wi(i)*wi(j)*pdf_ij*((sin(pi*(x(j)-c(i,j))/(T1_mode_1(i,j))))^2 + (sin(pi*(x(j)-c(i,j))/(T2_mode_1(i,j))))^2 + (sin(pi*(x(j)-c(i,j))/(T1_mode_2(i,j))))^2 + (sin(pi*(x(j)-c(i,j))/(T2_mode_2(i,j))))^2 + (sin(pi*(x(j)-c(i,j))/(T1_mode_3(i,j))))^2 + (sin(pi*(x(j)-c(i,j))/(T2_mode_3(i,j))))^2) + Int_eva;
        
    end
end

%The value of the scattering integral I
Int_eva=p1*q1*Int_eva;
end