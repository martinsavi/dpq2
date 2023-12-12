clear; clc;
a0 = [0.2,  0.2,1 , 0.05,2e-2, 45, 3.6e-6, 100]; % alfa a,alfa c, ecell, kk, rcell, ilim, Dglicerol
lb = [0.1 ,0.10,0 ,0.05, 2e-2, 45,  3.6e-6, 0];
ub = [0.5 ,0.50,2 ,2,    2e-2, 45,   3.6e-6, 100];

Vexp = [0.719347783
0.67935956
0.630218069
0.52992705
0.429229603
0.329989908
0.232060892
0.080883602
]'; % V
Iexp = [17.27941176 
18.97058824
21.17647059
23.45588235
27.20588235
31.32352941
35.66176471
42.86764706
]'; % A/m2
    v0 = [100 80 70 60 40 30 20 19 8.4 8.1 7.5 7.0 0.05 -0.2 0.4 6.0]'; %chute ini
    
options = optimset('TolFun',1e-20,'TolX',1e-20);
options.Algorithm = 'levenberg-marquardt';
% options.Algorithm = 'trust-region-reflective';
[x,fval] =  lsqnonlin(@(a)obj(a),a0,lb,ub,options);
Vteo = linspace(0.05,Vexp(1),30);

for i = 1:size(Vteo,2)
        options = optimset('Display','off');
        [sol] = fsolve(@(v) sistemaNL(x,v,Vteo(i)),v0,options);
        Iteo(i) = real(sol(16));
end

plot(Iexp,Vexp,'o',Iteo,Vteo,'r')
xlabel('I (A/m2)')
ylabel('Vcell (V)')
title(strcat('Curva desempenho - Fval = ',string(fval)))

function o = obj(a)
    
Vexp = [0.719347783
0.67935956
0.630218069
0.52992705
0.429229603
0.329989908
0.232060892
0.080883602
]'; % V
Iexp = [17.27941176 
18.97058824
21.17647059
23.45588235
27.20588235
31.32352941
35.66176471
42.86764706
]'; % A/m2
    Iteo = zeros(size(Vexp));
    v0 = [100 80 70 60 40 30 20 19 8.4 8.1 7.5 7.0 0.05 -0.2 0.4 6.0]'; %chute ini

    for i = 1:size(Vexp,1)
        options = optimset('Display','off');
        [sol] = fsolve(@(v) sistemaNL(a,v,Vexp(i)),v0,options);
        csvwrite('sols.csv',sol);
        Iteo(i) = sol(16);
    end
    c = [a+0.00]

   o = Iteo - Iexp;

      
end

function fv = sistemaNL(a,v,vcell)
    T = 303; %K
    R = 8.314; % m³.atm/K.mol
    P = 101325; %atm
    oxi0 = (0.21*P)/(R*T); % valor de entrada - mol/m³

    %D = 0.86^2.5*((T^1.75*5.8/10000)/(27.772*P))*3600/10000;
    D=a(7);
    K(1) = a(8);

    Iref=0.0064;
    Cref=1;
    
    ilim =a(6);
    v0 = 100 ;% mol/m³ - concentração inicial
    fx = 10 ;% assumido pelo artigo
    kk = a(4) ;% mol/m³.h 
    Km = a(8);
    K = [Km 1.0 0.8 0.8]'; % os 3 últimos valores foram assumidos pelo artigo (K2-4)
    alfaa = a(1); %0.5 % constante
    alfac = a(2); %0.44 % constante
    Kdec = 8.33*10^-4; % valor assumido pelo artigo
    q = 2.25*10^-5;  % vazão volumétrica - m³/h
    V = 1.596*10^-5; % m³ - volume
    F = 96485; % Constante de faraday - C/mol
    A = 1.2*10^-3; % área - m²
    YxA = 0.05; % rendimento 
    Ecell = a(3); %0.32; % volts
    Rcell = a(5);%4.9389*10^-5; % m³/S
    M=0.05; %quantidade de microorganismos no início (segundo artigo,0) mol/m³
    espce=0.000023; % assumido pelo artigo - m

  media_gli = (v(6)+v(7))/2;  
  fv = zeros(16,1);   
  fv(1) = v(1)- v0;
  fv(2) = v(2)- v(1);
  fv(3) = v(3)- v(2);
  fv(4) = v(4)- K(2)*v(3);
  fv(5) = v(5)- v(4);
  fv(6) = v(6)- K(3)*v(5);
  fv(7) = kk*exp((alfaa*v(14)*F)/R/T)*(media_gli)/(K(1)+media_gli)- D*(v(6)- v(7))/espce;
  fv(8) = v(8);
  fv(9) = oxi0 - v(9);
  fv(10)= v(10) - v(9);
  fv(11) = v(11)- K(4)*v(10);
  fv(12) = 3600*v(16)/4/F - (0.08*10^-4*3600)*(v(11)- v(12))/espce;
  fv(13) = v(13)-1;
  fv(14) = -vcell + Ecell -v(15)- v(16)*Rcell - v(14); 
  fv(15) = 3600*v(16)/14/F - kk*exp((alfaa*v(14)*F)/R/T)*media_gli/(K(1)+media_gli);
  fv(16) = v(16)-(Iref*(((v(11)+v(12))/2)/Cref)*exp(alfac*v(15)*F/R/T)); 
end



