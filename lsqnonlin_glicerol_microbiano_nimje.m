clear; clc;
a0 = [0.50  ,0.50  ,0.600,1,  4e-6 ,1];
lb = [0.10  ,0.10  ,0.600,1,4e-6 ,1 ];
ub = [0.90  ,0.90  ,0.600,1 , 4e-6,1 ];

Vexp = [0.523348018
0.522026432
0.500881057
0.453303965
0.413656388
0.397797357
0.383259912
0.362114537
0.340969163
0.321145374
0.292070485
0.27753304
0.27092511
0.264317181
0.255066079
0.251101322
0.249779736
]';
Iexp = [0
0.023006037
0.027319669
0.029476485
0.032352239
0.03738481
0.043855258
0.048887828
0.056077215
0.064704478
0.086991577
0.102808227
0.120781693
0.143068791
0.166793767
0.199146006
0.235092938
]';
    v0 = [22.26 15 13 12 11 10 5 2 8.4 8.1 7.5 7.0 0.05 -0.2 0.4 6.0]'; %chute ini
    
options = optimset('TolFun',1e-20,'TolX',1e-20);
options.Algorithm = 'levenberg-marquardt';
% options.Algorithm = 'trust-region-reflective';
[x,fval] =  lsqnonlin(@(a)obj(a),a0,lb,ub,options);
Vteo = linspace(0.25,0.6,30);

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
    
Vexp = [0.523348018
0.522026432
0.500881057
0.453303965
0.413656388
0.397797357
0.383259912
0.362114537
0.340969163
0.321145374
0.292070485
0.27753304
0.27092511
0.264317181
0.255066079
0.251101322
0.249779736
]';
Iexp = [0
0.023006037
0.027319669
0.029476485
0.032352239
0.03738481
0.043855258
0.048887828
0.056077215
0.064704478
0.086991577
0.102808227
0.120781693
0.143068791
0.166793767
0.199146006
0.235092938
]';
    Iteo = zeros(size(Vexp));
    v0 = [22.26 15 13 12 11 10 5 2 8.4 8.1 7.5 7.0 0.05 -0.2 0.4 6.0]'; %chute ini

    for i = 1:size(Vexp,1)
        options = optimset('Display','off');
        [sol] = fsolve(@(v) sistemaNL(a,v,Vexp(i)),v0,options);
        csvwrite('sols.csv',sol);
        Iteo(i) = sol(16);
    end
    c = [a+0.00]

    o = Iteo - Iexp;
%       o = sum((Iteo-Iexp).^2); %MMQ
      
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
    
    
    v0 = 1.56 ;% mol/m³ - concentração inicial
    fx = 10 ;% assumido pelo artigo
    kk = a(4) ;% mol/m³.h 
    K = [a(6) 0.8 0.8 0.8]'; % os 3 últimos valores foram assumidos pelo artigo (K2-4)
    alfaa = a(1); %0.5 % constante
    alfac = a(2); %0.44 % constante
    Kdec = 8.33*10^-4; % valor assumido pelo artigo
    q = 2.25*10^-7; % vazão volumétrica - m³/h
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
    fv(14) = -vcell + Ecell -v(15)- v(16)*Rcell - v(14); % primeira equação do frame 3
    fv(15) = 3600*v(16)/14/F - kk*exp((alfaa*v(14)*F)/R/T)*v(2)/(K(1)+v(2))*v(13); % terceira equação do frame 3
    fv(16) = v(16)-(Iref*(((v(11)+v(12))/2)/Cref)*exp(alfac*v(15)*F/R/T)); 
end



