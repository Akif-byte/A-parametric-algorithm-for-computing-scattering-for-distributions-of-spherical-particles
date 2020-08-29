%This script calculates |a_n|^2,|b_n|^2 using a naive
%discretization of the parametric space for different probability
%distributions p(m,x)=pdf_ij
clc
clear all
% Ranges for m1,m2 & x,y which covers the whole distribution
m1=[1.2:0.001:1.4];
x=[60:0.001:100];
m2=[1.51];
y=150;
%Bi-Modal Distribution Parameters
% mu1_x=75;mu2_x=80;sig1_x=10/3;sig2_x=10/3;
% mu1_m1=1.32;mu2_m1=1.28;sig1_m1=0.02;sig2_m1=0.02;
% p=0.5; %Weight for B-Modal distribution; f= p*Normal1 + (1-p)*Normal2
%Normal Distribution Parameters
% mu_x=80; mu_m1=1.3; sig_x=6.667; sig_m1=0.033;
Csca=0;

%Calculation of pdf_ij*(|a_n|^2+|b_n|^2) for the entire range of  the parametric
%space for the first 3 modes
for n=1:3
    nu=n+0.5;
    for j1=1:length(x)
        for j2=1:length(y)
            for k1=1:length(m1)
                for k2=1:length(m2)
                    z1=m1(k1)*x(j1);
                    z2=m2(k2)*x(j1);
                    z3=m2(k2)*y(j2);
                    sqy=sqrt(0.5*pi/y(j2));
                    sqz1=sqrt(0.5*pi/z1);
                    sqz2=sqrt(0.5*pi/z2);
                    sqz3=sqrt(0.5*pi/z3);
                    by=besselj(nu,y(j2))*sqy;
                    bz1=besselj(nu,z1)*sqz1;
                    bz2=besselj(nu,z2)*sqz2;
                    bz3=besselj(nu,z3)*sqz3;
                    yy=bessely(nu,y(j2))*sqy;
                    yz2=bessely(nu,z2)*sqz2;
                    yz3=bessely(nu,z3)*sqz3;
                    hy=by+1i*yy;
                    by_n_1=besselj(nu-1,y(j2))*sqy;
                    bz1_n_1=besselj(nu-1,z1)*sqz1;
                    bz2_n_1=besselj(nu-1,z2)*sqz2;
                    bz3_n_1=besselj(nu-1,z3)*sqz3;
                    yy_n_1=bessely(nu-1,y(j2))*sqy;
                    yz2_n_1=bessely(nu-1,z2)*sqz2;
                    yz3_n_1=bessely(nu-1,z3)*sqz3;
                    hy_n_1=by_n_1+1i*yy_n_1;
                    ay=y(j2)*by_n_1-n*by;
                    az1=z1*bz1_n_1-n*bz1;
                    az2=z2*bz2_n_1-n*bz2;
                    az3=z3*bz3_n_1-n*bz3;
                    cz2=z2*yz2_n_1-n*yz2;
                    cz3=z3*yz3_n_1-n*yz3;
                    ahy=y(j2)*hy_n_1-n*hy;
                    An=(m2(k2)*m2(k2)*bz2*az1-m1(k1)*m1(k1)*az2*bz1)/(m2(k2)*m2(k2)*yz2*az1-m1(k1)*m1(k1)*cz2*bz1);
                    Bn=(bz1*az2-bz2*az1)/(bz1*cz2-yz2*az1);
                    an=(by*(az3-An*cz3)-m2(k2)*m2(k2)*ay*(bz3-An*yz3))/(hy*(az3-An*cz3)-m2(k2)*m2(k2)*ahy*(bz3-An*yz3));
                    bn=(by*(az3-Bn*cz3)-ay*(bz3-Bn*yz3))/(hy*(az3-Bn*cz3)-ahy*(bz3-Bn*yz3));
                  
                    %Uniform Distribution
                    pdf_ij=1/((m1(end)-m1(1))*(x(end)-x(1)));
                   
                    %Normal Distribution
                    %pdf_ij=1/(2*pi*sig_x*sig_m1)*exp(-0.5*((((x(j1)-mu_x)/sig_x)^2) + ((m1(k1)-mu_m1)/sig_m1)^2));
                    
                    %Bimodal Distribution
%                     pdf1_ij=1/(2*pi*sig1_x*sig1_m1)*exp(-0.5*((((x(j1)-mu1_x)/sig1_x)^2) + ((m1(k1)-mu1_m1)/sig1_m1)^2)); 
%                     pdf2_ij=1/(2*pi*sig2_x*sig2_m1)*exp(-0.5*((((x(j1)-mu2_x)/sig2_x)^2) + ((m1(k1)-mu2_m1)/sig2_m1)^2));
%                     pdf_ij(j1,k1)=p*pdf1_ij+(1-p)*pdf2_ij;
                    Csca=pdf_ij*(2*n+1)*(abs(an)^2+abs(bn)^2)+Csca;
                end
            end
        end
    end
end

Csca=Csca*0.001^2;
