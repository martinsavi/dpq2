clear; clc;
a0=[0.50;0.44];

[x,fval] =  fmincon(@(a) obj(a),a0,[],[],[],[],[],[],@(a) cons(a));
% [x,fval] =  fmincon(@obj,x0,[],[],[],[],[],[],@cons);

function o = obj(a)
    Iexp = [.7/2,1.8/2,2.1/2,2.5/2,4.2/2,9/2,11/2,11.2/2]';
    Vexp = [.6/2,.53/2,.50/2,.49/2,.40/2,.20/2,.13/2,.10/2]';
    Iteo = zeros(size(Vexp));

    v0 = [1.25 1.20 1.15 1.0 0.9 0.8 0.7 0.66 8.4 8.1 7.5 7.0 0.05 -0.2 0.4 6.0]'; %chute ini
    
    
    for i = 1:size(Vexp,1)
        options = optimset('Display','off');
        [sol] = fsolve(@(v) sistemaNL(a,v,Vexp(i)),v0,options); 
        Iteo(i) = sol(16);
    end
    c = [Iteo+0.01; a+0.01]
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

    v0 = 1.56 ;% mol/m³ - concentração inicial
    fx = 10 ;% assumido pelo artigo
    kk = 0.207 ;% mol/m³.h 
    K = [0.592 0.8 0.8 0.8]'; % os 3 últimos valores foram assumidos pelo artigo (K2-4)
    alfaa = a(1); %0.5 % constante
    alfac = a(2); %0.44 % constante
    Kdec = 8.33*10^-4; % valor assumido pelo artigo
    q = 2.25*10^-5;  % vazão volumétrica - m³/h
    V = 1.596*10^-5; % m³ - volume
    F = 96485; % Constante de faraday - C/mol
    A = 1.2*10^-3; % área - m²
    YxA = 0.05; % rendimento 
    Ecell = 0.32; % volts
    Rcell = 4.9389*10^-5; % m³/S
    M=0.05; %quantidade de microorganismos no início (segundo artigo,0) mol/m³
    espce=0.000023; % assumido pelo artigo - m

    % v(1) = y(1)
    % v(2) = y(2)
    % v(3) = y(3)
    % v(4) = y(4)
    % v(5) = y(5)
    % v(6) = y(6)
    % v(7) = y(7) 
    % v(8) = y(8)
    % v(9) = oxi(1)
    % v(10) = oxi(2)
    % v(11) = oxi(3)
    % v(12) = oxi(4)
    % v(13) = microorganismo final
    % v(14) = etaA
    % v(15) = etaC
    % v(16)= Icell
    % YxA = rendimento
    % fx = wash out
    %fv = zeros(16,1);   % dimensionando fv
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
    fv(14) = -vcell + Ecell -v(15)- v(16)*Rcell - v(14)-R*T/F*log(6.45/(6.45-v(16))); % primeira equação do frame 3
    fv(15) = 3600*v(16)/20/F - kk*exp((alfaa*v(14)*F)/R/T)*v(2)/(K(1)+v(2))*v(13); % terceira equação do frame 3
    fv(16) = v(16)-(Iref*(((v(11)+v(12))/2)/Cref)*exp(alfac*v(15)*F/R/T)); 
end

function [c,ceq] = cons(a)

    T = 303; %K
    R = 8.314; % m³.atm/K.mol
    P = 101325; %atm
    oxi0 = (0.21*P)/(R*T); % valor de entrada - mol/m³

    %D = 0.86^2.5*((T^1.75*5.8/10000)/(27.772*P))*3600/10000;
    D=0.08*10^-4*3600;

    Iref=0.0064;
    Cref=1;

    v = [1.25 1.20 1.15 1.0 0.9 0.8 0.7 0.66 8.4 8.1 7.5 7.0 0.05 -0.2 0.4 6.0]'; %chute ini
    
    v0 = 1.56 ;% mol/m³ - concentração inicial
    fx = 10 ;% assumido pelo artigo
    kk = 0.207 ;% mol/m³.h 
    K = [0.592 0.8 0.8 0.8]'; % os 3 últimos valores foram assumidos pelo artigo (K2-4)
    alfaa = a(1); %0.5 % constante
    alfac = a(2); %0.44 % constante
    Kdec = 8.33*10^-4; % valor assumido pelo artigo
    q = 2.25*10^-5;  % vazão volumétrica - m³/h
    V = 1.596*10^-5; % m³ - volume
    F = 96485; % Constante de faraday - C/mol
    A = 1.2*10^-3; % área - m²
    YxA = 0.05; % rendimento 
    Ecell = 0.32; % volts
    Rcell = 4.9389*10^-5; % m³/S
    M=0.05; %quantidade de microorganismos no início (segundo artigo,0) mol/m³
    espce=0.000023; % assumido pelo artigo - m

    c=kk*exp((alfaa*v(14)*F)/R/T)*v(2)/(K(1)+v(2))*v(13)-(q/A)*(v(1)-v(2)); % Equação 3-13;
    ceq=[];
end
