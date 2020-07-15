%This function calculates the Average periods of sin^2(alpha) and
%sin^(beta) for a given mode and parameter c 
function [Avg_period_alpha,Avg_period_beta]=calc_period(c,n,m)
%Caluclation of sin^(alpha) and sin^(beta) along c=mx curves
nu=n+0.5;
x=c./m;
z=m.*x;
m2=m.*m;
sqx=sqrt(0.5*pi./x);
sqz=sqrt(0.5*pi./z);
bx=besselj(nu,x).*sqx;
bz=besselj(nu,z).*sqz;
yx=bessely(nu,x).*sqx;
hx=bx+1i*yx;
bx_n_1=besselj(nu-1,x).*sqx;
bz_n_1=besselj(nu-1,z).*sqz;
yx_n_1=bessely(nu-1,x).*sqx;
hx_n_1=bx_n_1+1i*yx_n_1;
ax=x.*bx_n_1-n*bx;
az=z.*bz_n_1-n*bz;
ahx=x.*hx_n_1-n*hx;
an=(m2.*bz.*ax-bx.*az)./(m2.*bz.*ahx-hx.*az);
bn=(bz.*ax-bx.*az)./(bz.*ahx-hx.*az);
tan_alpha_an = real(1i*an./(an-1));
tan_beta_bn = real(1i*bn./(bn-1));
alpha = atand(tan_alpha_an);
beta = atand(tan_beta_bn);      
sin_sq_alpha=sind(alpha).^2;
sin_sq_beta=sind(beta).^2;
count1=1;
count2=1;

%Finding the minima/maxima points of sin^(alpha) and sin^(beta) and
%calculating the mean of the difference between successive minima/maxima to
%get the average period 
if(sin_sq_alpha(length(sin_sq_alpha))<1e-6)
    x_min_alpha(count1)=x(length(x));
    count1=count1+1;
end
if(sin_sq_beta(length(sin_sq_beta))<1e-6)
    x_min_beta(count2)=x(length(x));
    count2=count2+1;
end
for i=2:length(x)-1
    if(sin_sq_alpha(i+1)>sin_sq_alpha(i) && sin_sq_alpha(i-1)>sin_sq_alpha(i) && sin_sq_alpha(i)<1e-2)
        x_min_alpha(count1)=x(i);
        count1=count1+1;
    end
    if(sin_sq_beta(i+1)>sin_sq_beta(i) && sin_sq_beta(i-1)>sin_sq_beta(i) && sin_sq_beta(i)<1e-2)
        x_min_beta(count2)=x(i);
        count2=count2+1;
    end
end
Avg_period_alpha=mean(diff(sort([x_min_alpha,c])));
Avg_period_beta=mean(diff(sort([x_min_beta,c])));
end