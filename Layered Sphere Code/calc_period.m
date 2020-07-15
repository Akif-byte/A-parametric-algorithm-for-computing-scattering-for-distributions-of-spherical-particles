%This function calculates the average period/phase for sin^2(alpha) and
%sin^2(beta) using the inputs as the parameter c1=m1*x, index of inner
%sphere m1, mode number n
function [Avg_period_alpha,Avg_period_beta,Phase_alpha,Phase_beta]=calc_period(c1,m1,n)
x=c1./m1;
m2=1.51;
y=150;
nu=n+0.5;
z1=m1.*x;
z2=m2.*x;
z3=m2.*y;
sqy=sqrt(0.5*pi./y);
sqz1=sqrt(0.5*pi./z1);
sqz2=sqrt(0.5*pi./z2);
sqz3=sqrt(0.5*pi./z3);
by=besselj(nu,y).*sqy;
bz1=besselj(nu,z1).*sqz1;
bz2=besselj(nu,z2).*sqz2;
bz3=besselj(nu,z3).*sqz3;
yy=bessely(nu,y).*sqy;
yz2=bessely(nu,z2).*sqz2;
yz3=bessely(nu,z3).*sqz3;
hy=by+1i*yy;
by_n_1=besselj(nu-1,y).*sqy;
bz1_n_1=besselj(nu-1,z1).*sqz1;
bz2_n_1=besselj(nu-1,z2).*sqz2;
bz3_n_1=besselj(nu-1,z3).*sqz3;
yy_n_1=bessely(nu-1,y).*sqy;
yz2_n_1=bessely(nu-1,z2).*sqz2;
yz3_n_1=bessely(nu-1,z3).*sqz3;
hy_n_1=by_n_1+1i*yy_n_1;
ay=y.*by_n_1-n*by;
az1=z1.*bz1_n_1-n*bz1;
az2=z2.*bz2_n_1-n*bz2;
az3=z3.*bz3_n_1-n*bz3;
cz2=z2.*yz2_n_1-n*yz2;
cz3=z3.*yz3_n_1-n*yz3;
ahy=y.*hy_n_1-n*hy;
An=(m2.*m2.*bz2.*az1-m1.*m1.*az2.*bz1)./(m2.*m2.*yz2.*az1-m1.*m1.*cz2.*bz1);
Bn=(bz1.*az2-bz2.*az1)./(bz1.*cz2-yz2.*az1);
an=(by.*(az3-An.*cz3)-m2.*m2.*ay.*(bz3-An.*yz3))./(hy.*(az3-An.*cz3)-m2.*m2.*ahy.*(bz3-An.*yz3));
bn=(by.*(az3-Bn.*cz3)-ay.*(bz3-Bn.*yz3))./(hy.*(az3-Bn.*cz3)-ahy.*(bz3-Bn.*yz3));
tan_alpha_an=real(1i*an./(an-1));
tan_beta_bn=real(1i*bn./(bn-1));
alpha=atand(tan_alpha_an);
beta=atand(tan_beta_bn);

sin_sq_alpha=sind(alpha).^2;
sin_sq_beta=sind(beta).^2;
count1=1;
count2=1;
count3=0;
count4=0;
%Calculation of average period of sin^2(alpha), sin^2(beta) by finding the
%minima/maxima and averaging 
for i=2:length(x)-1
        if(sin_sq_alpha(i+1)<sin_sq_alpha(i) && sin_sq_alpha(i-1)<sin_sq_alpha(i))
            x_max_alpha(count1)=x(i);
            count1=count1+1;
            count3=1;
        end
        
        if(sin_sq_beta(i+1)<sin_sq_beta(i) && sin_sq_beta(i-1)<sin_sq_beta(i))
            x_max_beta(count2)=x(i);
            count2=count2+1;
            count4=1;
        end        
end

count1=1;
count2=1;
    for i=2:length(x)-1
            if(sin_sq_alpha(i+1)>sin_sq_alpha(i) && sin_sq_alpha(i-1)>sin_sq_alpha(i))
                x_min_alpha(count1)=x(i);
                count1=count1+1;
            end
        if(sin_sq_beta(i+1)>sin_sq_beta(i) && sin_sq_beta(i-1)>sin_sq_beta(i))
            x_min_beta(count2)=x(i);
            count2=count2+1;
        end
    end
if(count3==1 || count1>1)
    if(length(x_max_alpha)>1 && length(x_min_alpha)>1)
        Avg_period_alpha=0.5*(mean(abs((diff((x_max_alpha))))) + mean(abs((diff(([x_min_alpha]))))));
    end

else
    Avg_period_alpha=2.1;
end

if(count4==1 || count2>1)
    if(length(x_max_beta)>1 && length(x_min_beta)>1)
        Avg_period_beta=0.5*(mean(abs((diff((x_max_beta))))) + mean(abs((diff(([x_min_beta]))))));
    end
else
    Avg_period_beta=2.1;
end
%Calculating phase
t1=0.4;
t2=-0.4;
Phase_alpha=-2*pi*x_min_alpha(end)/Avg_period_alpha - acos(t1) + 2*pi;
Phase_beta=-2*pi*x_min_beta(end)/Avg_period_beta - acos(t2) + 2*pi;
end