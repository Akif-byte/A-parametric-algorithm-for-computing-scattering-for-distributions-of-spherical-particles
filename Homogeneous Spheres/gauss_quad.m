%Gauss quadrature for evauluation of the integral I 
%a,b limits of outer integral
%c,d limts of inner integral
%p(m,x)=pdf_ij the probability density function
function [Int_eva]=gauss_quad(a,b,c,d,wi,xi)
%wi,xi are Weights and nodes for Gauss quadrature we have used the same
%number of weights and nodes for both parameters m and x

% Converting limits of integration to -1 to 1
p1=(b-a)/2;
p2=(b+a)/2;
q1=(d-c)/2;
q2=(d+c)/2;
x=p1*xi+p2;
m=q1*xi+q2; 

%Probabibility Distribution: Uniform Distribution
pdf_ij=1;

%Normal Distribution Parameters
% mu_x=15; mu_m=1.5; sig_x=1.67; sig_m=0.1;

%Bi-Modal Distribution Parameters
% mu1_x=13;mu2_x=17; mu1_m=1.4;mu2_m=1.6; sig1_x=1;sig2_x=1; sig1_m=0.06;sig2_m=0.06; 
% p=0.5; %Weight for B-Modal distribution; f= p*Normal1 + (1-p)*Normal2

Int_eva=0;

%Iterating over all the quadrature points in m, x and modes n
for i=1:length(m)
    for j=1:length(x)
        for n=1:3
%           Calculating the scattering coefficients analytically on
%           the nodes of Gauss-Legendre quadrature 
            nu=n+0.5;
            z=m(i)*x(j);
            m2=m(i)*m(i);
            sqx=sqrt(0.5*pi/x(j));
            sqz=sqrt(0.5*pi/z);
            bx=besselj(nu,x(j))*sqx;
            bz=besselj(nu,z)*sqz;
            yx=bessely(nu,x(j))*sqx;
            hx=bx+1i*yx;
            bx_n_1=besselj(nu-1,x(j))*sqx;
            bz_n_1=besselj(nu-1,z)*sqz;
            yx_n_1=bessely(nu-1,x(j))*sqx;
            hx_n_1=bx_n_1+1i*yx_n_1;
            ax=x(j)*bx_n_1-n*bx;
            az=z*bz_n_1-n*bz;
            ahx=x(j)*hx_n_1-n*hx;
            an=(m2*bz*ax-bx*az)/(m2*bz*ahx-hx*az);
            bn=(bz*ax-bx*az)/(bz*ahx-hx*az);
%           Calculating the integral I=Int_eva given by
%           pdf_ij*(|a_n|^2=sin^2(alpha)+|b_n|^2=sin^2(beta))

%           Normal distribution
%           pdf_ij=1/(2*pi*sig_x*sig_m)*exp(-0.5*((((x(j)-mu_x)/sig_x)^2) + ((m(i)-mu_m)/sig_m)^2));
            
            %Bi-Modal distribution
%           pdf1_ij=1/(2*pi*sig1_x*sig1_m)*exp(-0.5*((((x(j)-mu1_x)/sig1_x)^2) + ((m(i)-mu1_m)/sig1_m)^2)); 
%           pdf2_ij=1/(2*pi*sig2_x*sig2_m)*exp(-0.5*((((x(j)-mu2_x)/sig2_x)^2) + ((m(i)-mu2_m)/sig2_m)^2));
%           pdf_ij=p*pdf1_ij+(1-p)*pdf2_ij;
            
%           Multiplying by the weights of Gauss-Legendre quadrature 
            Int_eva=wi(i)*wi(j)*pdf_ij*((abs(an))^2 + (abs(bn))^2) + Int_eva;

        end
    end
end
%Value of the Integral I  
Int_eva=p1*q1*Int_eva;
end