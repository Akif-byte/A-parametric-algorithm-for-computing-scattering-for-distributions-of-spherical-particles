clc
clear all
%Input the patch size over which the integral is to be performed
%For uniform distribution it will be range of the interval in m and x
%For normal distribution m2-m1=6*sigma_m, x2-x1=6*sigma_x
m1=1.2;
m2=1.8;
x1=10;
x2=20;
%Normalized Scattering cross-section                                        
%Total scattering cross-section=2pi/k^2 * Csca; Where k is the wavenumber
%k= 2pi/lambda, lambda is the wavelength of incident light

%Loading Average periods for alpha, beta for different modes
% Alternatively the vectors Coeff_vec1, Coeff_vec2 can be loaded and the
% periods can be calucalted from them
load('Avg_period_alpha_mode_1')
load('Avg_period_beta_mode_1')
load('Avg_period_alpha_mode_2')
load('Avg_period_beta_mode_2')
load('Avg_period_alpha_mode_3')
load('Avg_period_beta_mode_3')

%Loading the precomuted weights and nodes for the gaussian quadrature. The
%weights are given by wi and the nodes are given by xi. We use the
%same number of nodes and weights for both m and x
load('GQ_x_i_20_w_i_20')

%Csca_gq computes the scattering integral using gaussian quadrature using
%the full bessel function evaluation of |a_n|^2, |b_n|^2

%Csca_algo computes the scattering integral using gaussian quadrature using
%the approximate trigonometric forms of |a_n|^2, |b_n|^2
% [Csca_gq]=gauss_quad(x1,x2,m1,m2,wi,xi);
[Csca_algo]=gauss_quad_appx_algo(x1,x2,m1,m2,Avg_period_alpha_mode_1,Avg_period_beta_mode_1,Avg_period_alpha_mode_2,Avg_period_beta_mode_2,Avg_period_alpha_mode_3,Avg_period_beta_mode_3,wi,xi);



