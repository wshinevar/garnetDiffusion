function [D]=calculateDiffusionMatrix(T,P,xalm,xpyr,xgrs,xsps)
%%
%this function calculates the diffusion matrix for garnet given the
%parameters from Borinski et al (2012) and Chu & Ague (2015) at fmq oxygen
%fugacity.

%The fo2 is calculated from fmq of oneill 1987a (huebner 1971 for
%pressure term) and cco= jakobsson & oskarsson 1994.

%T is temperature in kelvin
%P is pressure in bars
%xpyr... are the fraction of pyrope, almandine, grossular, and
%spessartine in the investigated garnet (they should add to 1)

%D is the multicomponent diffusion tensor [m^2/s]

%%

%calculate oxygen fugacity correction


R=8.31441;
fmq=(-587474+1584.7*T-203.3164*T*log(T)+0.09271*T^2+0.050*(P-1)/T)/R/T/log(10);
cco=4.325-21803/T+.171*(P-1)/T;
fo2correction=(fmq-cco)/6;


%calculate the mean ionic radius of the diffusive spot in the garnet
apyr=1.1456;
asps=1.1614;
agrs=1.1825;
aalm=1.1525;
a0=xalm*aalm+xgrs*agrs+xpyr*apyr+xsps*asps;



%all the following values are from Chu & Ague table 3 using equation 3
Dalm=10^(-7.98+457.2*(a0-aalm)/log(10)-(250.8e3+P*1.4)/R/T/log(10)+fo2correction);
Dpyr=10^(-8.28+340.3*(a0-aalm)/log(10)-(258e3+P*0.72)/R/T/log(10)+fo2correction);
Dsps=10^(-9.86+686.7*(a0-aalm)/log(10)-(212.9e3+P*0.95)/R/T/log(10)+fo2correction);
Dgrs=10^(-6.36+302.9*(a0-aalm)/log(10)-(299.3e3+P*1.77)/R/T/log(10)+fo2correction);


D=zeros(3);%almandine is 1, pyrope is 2, grossular 3, spessartine 4
%equation 1 in borinski et al 2012, this is the ideal case
D(1,1)=Dalm-Dalm*xalm/(Dalm*xalm+Dpyr*xpyr+Dgrs*xgrs+Dsps*xsps)*(Dalm-Dsps);
D(1,2)=-Dalm*xalm/(Dalm*xalm+Dpyr*xpyr+Dgrs*xgrs+Dsps*xsps)*(Dpyr-Dsps);
D(1,3)=-Dalm*xalm/(Dalm*xalm+Dpyr*xpyr+Dgrs*xgrs+Dsps*xsps)*(Dgrs-Dsps);

%pyrope
D(2,1)=-Dpyr*xpyr/(Dalm*xalm+Dpyr*xpyr+Dgrs*xgrs+Dsps*xsps)*(Dalm-Dsps);
D(2,2)=Dpyr-Dpyr*xpyr/(Dalm*xalm+Dpyr*xpyr+Dgrs*xgrs+Dsps*xsps)*(Dpyr-Dsps);
D(2,3)=-Dpyr*xpyr/(Dalm*xalm+Dpyr*xpyr+Dgrs*xgrs+Dsps*xsps)*(Dgrs-Dsps);

%grossular
D(3,1)=-Dgrs*xgrs/(Dalm*xalm+Dpyr*xpyr+Dgrs*xgrs+Dsps*xsps)*(Dalm-Dsps);
D(3,2)=-Dgrs*xgrs/(Dalm*xalm+Dpyr*xpyr+Dgrs*xgrs+Dsps*xsps)*(Dpyr-Dsps);
D(3,3)=Dgrs-Dgrs*xgrs/(Dalm*xalm+Dpyr*xpyr+Dgrs*xgrs+Dsps*xsps)*(Dgrs-Dsps);

end