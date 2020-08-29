%Gauss quadrature for evauluation of the integral I using the trigonometric approximation 
%a,b limits of outer integral
%c,d limts of inner integral
%p(m,x)=pdf_ij the probability density function
function [Csca]=gauss_quad_appx_algo(a1,b1,a2,b2,c1,d1,c2,d2,Avg_period_alpha_mode_1,Avg_period_beta_mode_1,Phase_alpha_mode_1,Phase_beta_mode_1,Avg_period_alpha_mode_2,Avg_period_beta_mode_2,Phase_alpha_mode_2,Phase_beta_mode_2,Avg_period_alpha_mode_3,Avg_period_beta_mode_3,Phase_alpha_mode_3,Phase_beta_mode_3,wi,xi)


%wi,xi are Weights and nodes for Gauss quadrature we have used the same
%number of weights and nodes for both parameters m1 and x
%x_i, m_i denotes the nodes for variables m and x respectively, i.e. roots
%of the legendre polynomial of degree 60
c1_i=wi;
c3_i=wi;
x1_i=xi;
m1_i=xi;

p1=(b1-a1)/2;
p2=(b1+a1)/2;
q1=(d1-c1)/2;
q2=(d1+c1)/2;

% Converting limits of intergration to -1 to 1
x=p1*x1_i+p2;
m1=q1*m1_i+q2;

%Uniform Distribution
pdf_ij=1/((b1-a1)*(d1-c1));

%Normal Distribution Parameters
% mu_x=80; mu_m1=1.3; sig_x=6.667; sig_m1=0.033;

%Bi-Modal Distribution Parameters
% mu1_x=75;mu2_x=80;sig1_x=10/3;sig2_x=10/3;
% mu1_m1=1.32;mu2_m1=1.28;sig1_m1=0.02;sig2_m1=0.02;
% p=0.5; %Weight for B-Modal distribution; f= p*Normal1 + (1-p)*Normal2

Int_eva=0;

%Scaling done to sin(Tx+phi) to fit it from 0 to 1 from -1 to 1
%t1 decides the degree of slantness
t1=0.4;
t2=-0.4;

%Maximum and minumum values of the function
%f_max,min=+-arctan(t/sqrt(1-t^2))/t
max_scale=2.05758;
min_scale=-1.02879;

%Calculating values of Average period, Phase for a given c=m1*x,
%Average_period_alpha/beta(1) corresponds to c=70,
%Average_period_alpha/beta(end) corresponds to c=170
c=m1*x';
tic
Index=floor((c-70)*1000+1);
period1_mode_1=Avg_period_alpha_mode_1(Index);
period2_mode_1=Avg_period_beta_mode_1(Index);
phase1_mode_1=Phase_alpha_mode_1(Index);
phase2_mode_1=Phase_beta_mode_1(Index);
period1_mode_2=Avg_period_alpha_mode_2(Index);
period2_mode_2=Avg_period_beta_mode_2(Index);
phase1_mode_2=Phase_alpha_mode_2(Index);
phase2_mode_2=Phase_beta_mode_2(Index);
period1_mode_3=Avg_period_alpha_mode_3(Index);
period2_mode_3=Avg_period_beta_mode_3(Index);
phase1_mode_3=Phase_alpha_mode_3(Index);
phase2_mode_3=Phase_beta_mode_3(Index);
toc
%gaussian quadrature over the entire parametric space in m and x but using
%the approximate trigonometric function forms instead of the bessel
%function form of |a_n|^2 & |b_n|^2
for i=1:length(m1)
    for j=1:length(x)
        %Calculating the periods and scaling to be done
        u1=2*pi*x(j)/period1_mode_1(i,j) + phase1_mode_1(i,j);
        u2=2*pi*x(j)/period2_mode_1(i,j) + phase2_mode_1(i,j);
        v1=2*pi*x(j)/period1_mode_2(i,j) + phase1_mode_2(i,j);
        v2=2*pi*x(j)/period2_mode_2(i,j) + phase2_mode_2(i,j);
        w1=2*pi*x(j)/period1_mode_3(i,j) + phase1_mode_3(i,j);
        w2=2*pi*x(j)/period2_mode_3(i,j) + phase2_mode_3(i,j);
        
        sin_sq_alpha_appx_mode_1=(atan((t1*sin(u1))/(1-t1*cos(u1))))/t1;
        sin_sq_beta_appx_mode_1=(atan((t2*sin(u2))/(1-t2*cos(u2))))/t2;
        sin_sq_alpha_appx_mode_2=(atan((t2*sin(v1))/(1-t2*cos(v1))))/t2;
        sin_sq_beta_appx_mode_2=(atan((t1*sin(v2))/(1-t1*cos(v2))))/t1;
        sin_sq_alpha_appx_mode_3=(atan((t1*sin(w1))/(1-t1*cos(w1))))/t1;
        sin_sq_beta_appx_mode_3=(atan((t2*sin(w2))/(1-t2*cos(w2))))/t2;

        sin_sq_alpha_appx_mode_1=sin_sq_alpha_appx_mode_1-min_scale;
        sin_sq_alpha_appx_mode_1=sin_sq_alpha_appx_mode_1/max_scale;
        sin_sq_beta_appx_mode_1=sin_sq_beta_appx_mode_1-min_scale;
        sin_sq_beta_appx_mode_1=sin_sq_beta_appx_mode_1/max_scale;

        sin_sq_alpha_appx_mode_2=sin_sq_alpha_appx_mode_2-min_scale;
        sin_sq_alpha_appx_mode_2=sin_sq_alpha_appx_mode_2/max_scale;
        sin_sq_beta_appx_mode_2=sin_sq_beta_appx_mode_2-min_scale;
        sin_sq_beta_appx_mode_2=sin_sq_beta_appx_mode_2/max_scale;
        sin_sq_alpha_appx_mode_3=sin_sq_alpha_appx_mode_3-min_scale;
        sin_sq_alpha_appx_mode_3=sin_sq_alpha_appx_mode_3/max_scale;
        sin_sq_beta_appx_mode_3=sin_sq_beta_appx_mode_3-min_scale;
        sin_sq_beta_appx_mode_3=sin_sq_beta_appx_mode_3/max_scale;
 
        %Normal Distribution
%         pdf_ij=1/(2*pi*sig_x*sig_m1)*exp(-0.5*((((x(j)-mu_x)/sig_x)^2) + ((m1(i)-mu_m1)/sig_m1)^2));

%       Bimodal Distribution
%         pdf1_ij=1/(2*pi*sig1_x*sig1_m1)*exp(-0.5*((((x(j)-mu1_x)/sig1_x)^2) + ((m1(i)-mu1_m1)/sig1_m1)^2)); 
%         pdf2_ij=1/(2*pi*sig2_x*sig2_m1)*exp(-0.5*((((x(j)-mu2_x)/sig2_x)^2) + ((m1(i)-mu2_m1)/sig2_m1)^2));
%         pdf_ij=p*pdf1_ij+(1-p)*pdf2_ij;
        
        %Multiplying by the weights for the Gauss-Legendre quadrature
        Int_eva=c1_i(i)*c3_i(j)*pdf_ij*(3*sin_sq_alpha_appx_mode_1 + 3*sin_sq_beta_appx_mode_1 + 5*sin_sq_alpha_appx_mode_2 + 5*sin_sq_beta_appx_mode_2 + 7*sin_sq_alpha_appx_mode_3 + 7*sin_sq_beta_appx_mode_3) + Int_eva; 
    end
end
toc
%The value of the scattering integral I
Csca=p1*q1*Int_eva;
end


