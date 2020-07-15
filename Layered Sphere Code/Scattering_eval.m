clc
clear all
%Input the patch size over which the integral is to be performed
%For uniform distribution it will be range of the interval in m1 and x
m1=1.2;
m2=1.4;
m3=1.51;
m4=1.51;
x1=60;
x2=100;
x3=150;
x4=150;
%Normalized Scattering cross-section                                        
%Total scattering cross-section=2pi/k^2 * Csca; Where k is the wavenumber
%k= 2pi/lambda, lambda is the wavelength of incident light

%Loading Average periods & phase for alpha, beta for different modes
load('Avg_period_alpha_mode_1')
load('Avg_period_beta_mode_1')
load('Phase_alpha_mode_1')
load('Phase_beta_mode_1')
load('Avg_period_alpha_mode_2')
load('Avg_period_beta_mode_2')
load('Phase_alpha_mode_2')
load('Phase_beta_mode_2')
load('Avg_period_alpha_mode_3')
load('Avg_period_beta_mode_3')
load('Phase_alpha_mode_3')
load('Phase_beta_mode_3')

%Loading the precomuted weights and nodes for the gaussian quadrature. The
%weights are given by wi and the nodes are given by xi. We use the
%same number of nodes and weights for both m and x
load('GQ_x_i_15_w_i_15')

%Csca_gq computes the scattering integral using gaussian quadrature using
%the full bessel function evaluation of |a_n|^2, |b_n|^2

%Csca_algo computes the scattering integral using gaussian quadrature using
%the approximate trigonometric forms of |a_n|^2, |b_n|^2
[Csca_gq]=gauss_quad_2D(x1,x2,x3,x4,m1,m2,m3,m4,wi,xi);    
[Csca_alg]=gauss_quad_appx_algo(x1,x2,x3,x4,m1,m2,m3,m4,Avg_period_alpha_mode_1,Avg_period_beta_mode_1,Phase_alpha_mode_1,Phase_beta_mode_1,Avg_period_alpha_mode_2,Avg_period_beta_mode_2,Phase_alpha_mode_2,Phase_beta_mode_2,Avg_period_alpha_mode_3,Avg_period_beta_mode_3,Phase_alpha_mode_3,Phase_beta_mode_3,wi,xi);












