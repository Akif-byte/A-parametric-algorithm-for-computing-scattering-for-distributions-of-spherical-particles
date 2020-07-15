%Gauss quadrature for 2D integration a,b limits of outer integral
%c,d limts of inner integral
%The integral is evaulated over m1 and x parametric space
%p(m,x)=pdf_ij the probability density function
function [Int_eva]=gauss_quad_2D(a1,b1,a2,b2,c1,d1,c2,d2,wi,xi)
%Weights and nodes for Gauss quadrature are wi,xi

c1_i=wi;
c3_i=wi;
x1_i=xi;
m1_i=xi;

p1=(b1-a1)/2;
p2=(b1+a1)/2;
q1=(d1-c1)/2;
q2=(d1+c1)/2;

%Converting limits of intergration to -1 to 1
x=p1*x1_i+p2;
m1=q1*m1_i+q2;
y=a2;   
m2=c2;
%Uniform Distribution
pdf_ij=1;

%Normal Distribution Parameters
% mu_x=80; mu_m1=1.3; sig_x=6.667; sig_m1=0.033;

%Bi-Modal Distribution Parameters
% mu1_x=75;mu2_x=80;sig1_x=10/3;sig2_x=10/3;
% mu1_m1=1.32;mu2_m1=1.28;sig1_m1=0.02;sig2_m1=0.02;
% p=0.5; %Weight for B-Modal distribution; f= p*Normal1 + (1-p)*Normal2

Int_eva=0;
%Iterating over all the quadrature points in m, x and modes n
for j1=1:length(x)
    for k1=1:length(m1)
        for n=1:3 
            nu=n+0.5;
            z1=m1(k1)*x(j1);
            z2=m2*x(j1);
            z3=m2*y;
            sqy=sqrt(0.5*pi/y);
            sqz1=sqrt(0.5*pi/z1);
            sqz2=sqrt(0.5*pi/z2);
            sqz3=sqrt(0.5*pi/z3);
            by=besselj(nu,y)*sqy;
            bz1=besselj(nu,z1)*sqz1;
            bz2=besselj(nu,z2)*sqz2;
            bz3=besselj(nu,z3)*sqz3;
            yy=bessely(nu,y)*sqy;
            yz2=bessely(nu,z2)*sqz2;
            yz3=bessely(nu,z3)*sqz3;
            hy=by+1i*yy;
            by_n_1=besselj(nu-1,y)*sqy;
            bz1_n_1=besselj(nu-1,z1)*sqz1;
            bz2_n_1=besselj(nu-1,z2)*sqz2;
            bz3_n_1=besselj(nu-1,z3)*sqz3;
            yy_n_1=bessely(nu-1,y)*sqy;
            yz2_n_1=bessely(nu-1,z2)*sqz2;
            yz3_n_1=bessely(nu-1,z3)*sqz3;
            hy_n_1=by_n_1+1i*yy_n_1;
            ay=y*by_n_1-n*by;
            az1=z1*bz1_n_1-n*bz1;
            az2=z2*bz2_n_1-n*bz2;
            az3=z3*bz3_n_1-n*bz3;
            cz2=z2*yz2_n_1-n*yz2;
            cz3=z3*yz3_n_1-n*yz3;
            ahy=y*hy_n_1-n*hy;
            An=(m2*m2*bz2*az1-m1(k1)*m1(k1)*az2*bz1)/(m2*m2*yz2*az1-m1(k1)*m1(k1)*cz2*bz1);
            Bn=(bz1*az2-bz2*az1)/(bz1*cz2-yz2*az1);
            an=(by*(az3-An*cz3)-m2*m2*ay*(bz3-An*yz3))/(hy*(az3-An*cz3)-m2*m2*ahy*(bz3-An*yz3));
            bn=(by*(az3-Bn*cz3)-ay*(bz3-Bn*yz3))/(hy*(az3-Bn*cz3)-ahy*(bz3-Bn*yz3));

            %Normal Distribution
%             pdf_ij=1/(2*pi*sig_x*sig_m1)*exp(-0.5*((((x(j1)-mu_x)/sig_x)^2) + ((m1(k1)-mu_m1)/sig_m1)^2));
            
            %Bimodal Distribution
            %pdf1_ij=1/(2*pi*sig1_x*sig1_m1)*exp(-0.5*((((x(j1)-mu1_x)/sig1_x)^2) + ((m1(k1)-mu1_m1)/sig1_m1)^2));
            %pdf2_ij=1/(2*pi*sig2_x*sig2_m1)*exp(-0.5*((((x(j1)-mu2_x)/sig2_x)^2) + ((m1(k1)-mu2_m1)/sig2_m1)^2));
            %pdf_ij=p*pdf1_ij+(1-p)*pdf2_ij;
            
            %Multiplying by the weights of Gauss-Legendre quadrature 
            Int_eva=c1_i(j1)*c3_i(k1)*pdf_ij*(abs(an)^2+abs(bn)^2) + Int_eva; 
        end 
    end
end
%Value of the Integral I  
Int_eva=p1*q1*Int_eva;
end