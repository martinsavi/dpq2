clear; clc;
a0 = [0.4575    0.1390    0.6661    0.1563    0.4000  319.9777  28.56   2.55e-5  ];
lb = [0.10     ,0.10     ,0.500    ,0.1      , 4e-2 ,100       , 1 , 1e-7  ];
ub = [0.9       ,0.2     ,0.700    ,0.5      , 4e-0 ,1000     , 100,  1e-3 ];

Vexp = [0.060048426
0.102663438
0.226634383
0.282808717
0.315738499
0.414527845
0.528813559
0.592736077
0.652784504
0.666343826
0.676029056

]';
Iexp = [0.38424015
0.354221388
0.267166979
0.225891182
0.206378987
0.141838649
0.074296435
0.042026266
0.015009381
0.009005629
0.005253283

]';
    v0 = [.125 .120 .115 .10 .09 .08 .07 .066 8.4 8.1 7.5 7.0 0.05 -0.2 0.4 0.6]'; %chute ini

%options = optimoptions('simulannealbnd','TimeLimit',60*5);
%[x,fval] = simulannealbnd(@(a)obj(a),a0,lb,ub,options);    

options = optimoptions('ga','ConstraintTolerance',1e-8,'PlotFcn', @gaplotbestf);
%options = optimset('TolFun',1e-20,'TolX',1e-20);
[x,fval] = ga(@(a)obj(a),8,[],[],[],[],lb,ub,[],options);

%options.Algorithm = 'levenberg-marquardt';
%[x,fval] =  lsqnonlin(@(a)obj(a),a0,lb,ub,options);

%options = optimset('TolFun',1e-20,'TolX',1e-20);
Vteo = linspace(Vexp(1),1,30);

for i = 1:size(Vteo,2)
        options = optimset('Display','off');
        [sol] = fsolve(@(v) sistemaNL(x,v,Vteo(i)),v0,options);
        Iteo(i) = real(sol(16));
end

plot(Iexp,Vexp,'o',Iteo,Vteo)
xlabel('I (A/m2)')
ylabel('Vcell (V)')
title(strcat('Curva desempenho - Fobj = ',string(fval)))

function o = obj(a)
    
    Vexp = [0.060048426
0.102663438
0.226634383
0.282808717
0.315738499
0.414527845
0.528813559
0.592736077
0.652784504
0.666343826
0.676029056
0.695399516
]';
Iexp = [0.38424015
0.354221388
0.267166979
0.225891182
0.206378987
0.141838649
0.074296435
0.042026266
0.015009381
0.009005629
0.005253283
0.000750469
]';
    Iteo = zeros(size(Vexp));
    v0 = [.125 .120 .115 .10 .09 .08 .07 .066 8.4 8.1 7.5 7.0 0.05 -0.2 0.4 0.6]'; %chute ini

    for i = 1:size(Vexp,1)
        options = optimset('Display','off');
        [sol] = fsolve(@(v) sistemaNL(a,v,Vexp(i)),v0,options);
        csvwrite('sols.csv',sol);
        Iteo(i) = sol(16);
    end
    c = [a+0.00]

    %o = Iteo - Iexp;
    o = sum((Iteo-Iexp).^2); %MMQ
    %o = sum(abs(Iteo-Iexp));
 
      
end

function fv = sistemaNL(a,v,vcell)
    T = 35.5+273.15; %K
    R = 8.314; % m³.atm/K.mol
    P = 101325; %atm
    oxi0 = (0.21*P)/(R*T); % valor de entrada - mol/m³

    %D = 0.86^2.5*((T^1.75*5.8/10000)/(27.772*P))*3600/10000;
    D=0.08*10^-4*3600;

    Iref=0.0064;
    Cref=1;
    
    
    v0 = a(7) ;% mol/m³ - concentração inicial
    fx = 10 ;% assumido pelo artigo
    kk = a(4) ;% mol/m³.h 
    K = [a(6) 0.8 0.8 0.8]'; % os 3 últimos valores foram assumidos pelo artigo (K2-4)
    alfaa = a(1); %0.5 % constante
    alfac = a(2); %0.44 % constante
    Kdec = 8.33*10^-4; % valor assumido pelo artigo
    q = a(8); % vazão volumétrica - m³/h
    V = 1.596*10^-5; % m³ - volume
    F = 96485; % Constante de faraday - C/mol
    A = 1.2*10^-3; % área - m²
    YxA = 0.05; % rendimento 
    Ecell = a(3); %0.32; % volts
    Rcell = a(5);%4.9389*10^-5; % m³/S
    M=0.05; %quantidade de microorganismos no início (segundo artigo,0) mol/m³
    espce=0.000023; % assumido pelo artigo - m

    fv = zeros(16,1);
    fv(1) = v(1)- v0;
    fv(2) = kk*exp((alfaa*v(14)*F)/R/T)*v(2)/(K(1)+v(2))*v(13)-(q/A)*(v(1)-v(2)); % Equação 3-13
    fv(3) = v(3)- v(2);
    fv(4) = v(4)- K(2)*v(3);
    fv(5) = v(5)- v(4);
    fv(6) = v(6)- K(3)*v(5);
    fv(7) = v(7)- v(6);
    fv(8) = v(8);
    fv(9) = oxi0 - v(9);
    fv(10)= v(10) - v(9);
    fv(11) = v(11)- K(4)*v(10);
    fv(12) = 3600*v(16)/4/F - D*(v(11)- v(12))/espce;
    % fv(13) = ((V*Kdec*v(13))/(A*YxA)) + ((q)/(A*fx*YxA)*(M - v(13))); % equação 4 do artigo igualada a zero
    fv(13) = kk*exp((alfaa*v(14)*F)/R/T)*v(2)/(K(1)+v(2))*v(13)-((V*Kdec*v(13))/(A*YxA)) + ((q)/(A*fx*YxA)*(M - v(13)));
    fv(14) = -vcell + Ecell - v(16)*Rcell-v(14)-v(15) ; % primeira equação do frame 3
    fv(15) = 3600*v(16)/14/F - kk*exp((alfaa*v(14)*F)/R/T)*v(2)/(K(1)+v(2))*v(13); % terceira equação do frame 3
    fv(16) = v(16)-(Iref*(((v(11)+v(12))/2)/Cref)*exp(alfac*v(15)*F/R/T)); 
end



