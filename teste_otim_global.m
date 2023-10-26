clear; clc;
a0 = [0.6000,0.4880,0.3453,0.1527,0.0005,4.2000];
lb = [0.45  ,0.30  ,0.20  ,0.150 ,4e-5   ,4.2];
ub = [0.80  ,0.70  ,0.50  ,0.550,5e-4,4.4  ];
 
melhor_ate_agora=[0.5592,    0.4500,    0.3565,    0.1560,    4e-5,    4.2130];
%   Iexp = [1.03759418890913	1.93984923218668	4.2105267859144	5.50375960737219	6.48120319236476	7.39849617344079	8.25563855060027	10.2857136886617	11.0676692450222	12.0902260697715
% ]';
%     Vexp = [0.66286288281833	0.612695560680013	0.51870282474702	0.517403067839343	0.48963963129269	0.462366184933421	0.423861324718456	0.326127045742691	0.279645383820309	0.229047044443503
% ]';  dados Borole-Fig4-Furfural - fica ruim

     Iexp = [0.4358974716599329	0.7179486905462856	1.7777778884712203	2.034188146739269	2.2222224013443386	2.585470267224074	2.8632480470144603	3.431623923883132	3.7905986767094046	3.9829063704104404
]';
    Vexp = [0.2933823488265998	0.26198529829719064	0.20330882067808234	0.18786764442405074	0.17654411124420005	0.15904412327031353	0.14463236269344493	0.11272059450821882	0.09727941825418726	0.07772059892597478
]'; % dados semi 2 DPQ 2 - Furfural

    v0 = [1.25 1.20 1.15 1.0 0.9 0.8 0.7 0.66 8.4 8.1 7.5 7.0 0.05 -0.2 0.4 6.0]'; %chute ini
    
% options = optimoptions('ga','TimeLimit',60*10);
% [x,fval] = ga(@(a) obj(a,Iexp,Vexp,v0),size(a0,2),[],[],[],[],lb,ub,[],options);

options = optimoptions('simulannealbnd','TimeLimit',60*2);
[x,fval] = simulannealbnd(@(a) obj(a,Iexp,Vexp,v0),a0,lb,ub,options);
Vteo = linspace(0.05,Vexp(1),30);

for i = 1:size(Vteo,2)
        options = optimset('Display','off');
        [sol] = fsolve(@(v) sistemaNL(x,v,Vteo(i)),v0,options);
        Iteo(i) = real(sol(16));
end

plot(Iexp,Vexp,'o',Iteo,Vteo)
xlabel('I (A/m2)')
ylabel('Vcell (V)')
title(strcat('Curva desempenho - Fobj = ',string(fval)))

function o = obj(a,Iexp,Vexp,v0)
    Iteo = zeros(size(Vexp));
    
    for i = 1:size(Vexp,1)
        options = optimset('Display','off');
        [sol] = fsolve(@(v) sistemaNL(a,v,Vexp(i)),v0,options);
        csvwrite('sols.csv',sol)
        Iteo(i) = sol(16);
    end
    c = [a+0.00001]
%       o = (1-a(1))^2 +100*(a(2)-a(1)^2)^2; Rosenbrock funciona ok
      o = sum((Iteo-Iexp).^2); %MMQ
      
end

function fv = sistemaNL(a,v,vcell)
    T = 303; %K
    R = 8.314; % m³.atm/K.mol
    P = 101325; %atm
    oxi0 = (0.21*P)/(R*T); % valor de entrada - mol/m³

    %D = 0.86^2.5*((T^1.75*5.8/10000)/(27.772*P))*3600/10000;
    D=0.08*10^-4*3600;

    Iref=0.0064;
    Cref=1;
    
    ilim =a(6);
    v0 = 1.56 ;% mol/m³ - concentração inicial
    fx = 10 ;% assumido pelo artigo
    kk = a(4) ;% mol/m³.h 
    K = [0.592 0.8 0.8 0.8]'; % os 3 últimos valores foram assumidos pelo artigo (K2-4)
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
    fv(14) = -vcell + Ecell -v(15)- v(16)*Rcell - v(14)-R*T/F*log(ilim/sqrt((ilim-v(16))^2)); % primeira equação do frame 3
    fv(15) = 3600*v(16)/20/F - kk*exp((alfaa*v(14)*F)/R/T)*v(2)/(K(1)+v(2))*v(13); % terceira equação do frame 3
    fv(16) = v(16)-(Iref*(((v(11)+v(12))/2)/Cref)*exp(alfac*v(15)*F/R/T)); 
end





