function [J, Dg,  ng] = JNNF_HET_Barahona2018(zeta, T, Swx)

%Heterogeneous nucleation rate in immersion mode following Barahona, ACP, 2018.
%Used in the final paper.
%input
% Swx is the water activity
% zeta is the templating factor
% T is the temperature in K 
%output
% J (1/(s*m2)), DG (work of nucleation J/K), 
% ng : germ size
% mode: 

%Constants
Tom=273.15;
R=8.314; %J/molK
Mw=18e-3; %Kg/mol
kB=1.3806e-23; %Boltzman constant J/K
Na=6.023e23; % Avogadro's constant molec/mol
E=892; %Smith and Kay, 1999'
T0=118; %Smith and Kay, 1999'
D0=3.06e-7; %Smith and Kay, 1999'
Tc05 =  211.473; %B18
z = 6; %Holten et al. 2013
Gw=1.46;  %from Spaepen 1975, assuming hcp crystals (Barahona 2014)
ss=1.105; %B14
nt=16; %B15


%Physical properties
aweq=awice_koop2009X(T);
Tr=(T-Tom)/Tom;
rom=IceDenX(Tom); % Kg/m3
vw1=1-(0.05294*Tr)-(0.05637*Tr*Tr)-(0.002913*Tr*Tr*Tr);
vw=Mw/(Na*rom*vw1); %m3 Zobrist et al 2007
vw0=Mw/(Na*rom); %m3 Zobrist et al 2007
S0=  (-7.7481e-5*T*T+5.5160e-2*T-6.6716)*R*vw/vw0; % configurational entropy (Barahona, 2018) 
ATS = E/(T-T0);
Dinfty=D0*exp(-ATS);
Dhm=DHX(T);


% Effective aw
DSmix =(zeta*log(zeta+1e-6) + (1-zeta)*log(1-zeta+ 1e-6))/z; %we add 1e-6 for numerical stability 
Aw  = 2*Tc05/(z*T);
Lmix=  Aw*(1-zeta)*zeta + DSmix;
aux=1/(1-zeta);
Sw=(Swx/(aweq^zeta))^aux;
Sw=Sw*exp(-Lmix/(1-zeta));

% Work of nucleation and germ size
Dg1=log(Sw/aweq) + log(Sw);
Dgb= kB*T*(Dg1);
Dgh=(Dhm-Gw*R*T*log(Sw))*ss*Gw/(Na);
ng0= (8/27)*((Dgh/Dgb)^3);
ng0=max(ng0, 0);
ng=ng0+2; %THis is to account for dissipated work
Dg=0.5*Dgb*ng;

% Interfacial diffusion
DSfideal = R*log(Swx/aweq);
Wdiss=(1-zeta)*kB*T*log(Swx/aweq);
f = exp(nt*Wdiss/(kB*T));
Sc = S0*(1-zeta) + zeta*DSfideal; %
g_fh =(S0/Sc)-1;
Ds= Dinfty*exp(-ATS*g_fh)*(1/(1+f));

% Nucleation rate
Zel=sqrt(Dg/(3*pi*kB*T*ng*ng));
d0=(6*vw/pi)^(1/3);
a0=0.25*pi*(d0^2);
n1=1/a0;
Ag=Gw*ss*(ng^(2/3))*a0;
fstar= Ag*Ds/(vw*d0);
Jprime =  Zel*fstar*n1;
J=Jprime*exp(-Dg/(kB*T));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function [aw] = awice_koop2009X(Te)
        % From Koop et.al. (2000), awice(T) (pure water) T in K
        imax=length(Te);
        aw_=Te*0;
        c1 = 210368;
        c2 = 131.438;
        c3 = -3323730;
        c4 = -41729.1;
        c5 = 8.31441;
        for i=1:imax    
            Tx=Te(i);
            aux1 = exp(9.550426 -(5723.265/Tx) + 3.53068*log(Tx) - 0.00728332*Tx);
            aux2 = 54.842763 -(6763.22/Tx) -4.210*log(Tx) + 0.000367*Tx; 
            aux3 = tanh(0.0415*(Tx-218.8))*(53.878 - (1331.22/Tx) - 9.44523*log(Tx) + 0.014025*Tx);
            aux4=exp(aux2+aux3);
            aw_(i) = aux1/aux4;
        end
        aw=aw_;
    end 
    
    function  [den]=IceDenX(Tx)
        %Ice density down to 93K in Kg/m3 from PK97.
        Tc = Tx - 273;
        a=zeros(3);
        a(1) = 0.9167;
        a(2) = -1.75e-4;
        a(3) = -5.0e-7;
        aux = 0;
        for i=1:3  
            aux= aux + (a(i) * Tc ^ (i - 1));
        end
        den=aux*1e3;
    end 

    function Dh = DHX(Te)
    %From Johari(1993)
        %J/mol
        Dh=Te*0.0;
        imax=length(Te);
        for i=1:imax
          Tx=Te(i);  
          Dh(i)=7.50856E-7*Tx^5 - 8.40025E-4*Tx^4 + 3.67171E-1*Tx^3 -  7.81467E1*Tx^2 ...
              + 8.11702E3*Tx - 3.29032E5; %Fit used in Barahona 2014, 2015..
        end 
    end 
end 
