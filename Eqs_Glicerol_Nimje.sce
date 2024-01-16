// DPQ2 - Professor Ruy de Sousa Júnior
// Leonardo Pacheco RA: 773350
// Marcus Vinicius Sgobi RA: 769819
// Maria Clara Sertori RA: 769750
// Paulo Eduardo da Silva Strozzi RA: 769926
// Stephanie Aline Souza Rodrigues RA: 770711
//for i =1:10; clf(i);end
clf(3)
function fv = funcv(v,v0,kk,K,alfaa,alfac,q,V,F,A,YxA,fx,espce,Kdec,Rcell,vcell,Ecell,M,D,T,ilim)
// Calcula o vetor de funções f(v), onde f(v[solução])=0 

R = 8.314 // m³.atm/K.mol
P = 101325 //atm
oxi0 = (0.21*P)/(R*T) // valor de entrada - mol/m³

//D = 0.86^2.5*((T^1.75*5.8/10000)/(27.772*P))*3600/10000;
//D=0.08*10^-4*3600

Iref=0.0064
Cref=1
//Iref = (4.222*10^-2)*exp((73200/R)*(1/353 - 1/T));
//Cref = (0.21*P)/(R*T);
  // v(1) = y(1)
  // v(2) = y(2)
  // v(3) = y(3)
  // v(4) = y(4)
  // v(5) = y(5)
  // v(6) = y(6)
  // v(7) = y(7) 
  // v(8) = y(8)
  // v(9) = oxi(1)
  // v(10) = oxi(2)
  // v(11) = oxi(3)
  // v(12) = oxi(4)
  // v(13) = microorganismo final
  // v(14) = etaA
  // v(15) = etaC
  // v(16)= Icell
  // YxA = rendimento
  // fx = wash out
  fv = zeros(16,1);   // dimensionando fv
  fv(1) = v(1)- v0;
  fv(2) = kk*exp((alfaa*v(14)*F)/R/T)*v(2)/(K(1)+v(2))*v(13)-(q/A)*(v(1)-v(2)); // Equação 3-13
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
  fv(13) = kk*exp((alfaa*v(14)*F)/R/T)*v(2)/(K(1)+v(2))*v(13)-((V*Kdec*v(13))/(A*YxA)) + ((q)/(A*fx*YxA)*(M - v(13)));
    fv(14) = -vcell + Ecell -v(15)- v(16)*Rcell - v(14)-R*T/F*log(ilim/(ilim-v(16))); // primeira equação do frame 3
    fv(15) = 3600*v(16)/14/F - kk*exp((alfaa*v(14)*F)/R/T)*v(2)/(K(1)+v(2))*v(13); // terceira equação do frame 3
    fv(16) = v(16)-(Iref*(((v(11)+v(12))/2)/Cref)*exp(alfac*v(15)*F/R/T)); 
endfunction

// Célula a combustivel microbiana
//clear // limpa memória do scilab
clc  // limpa a tela da janela de comando
Iexp = [0.029476485
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
0.268883055
0.290451214
0.311300435
]';
    Vexp = [0.453303965
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
0.236563877
0.216740088
0.188986784
]';

// Parâmetros do processo (sistema)
N = 16
T = 303//K
v0 = 1.56 ;// mol/m³ - concentração inicial
fx = 10 // assumido pelo artigo
v = [.125 .120 .115 .10 .09 .08 .07 .066 8.4 8.1 7.5 7.0 0.05 -0.2 0.4 0.6]'; // mol/m³
kk = 0.0000207/2 ;// mol/m³.h 
K = [0.592 0.8 0.8 0.8]'; // os 3 últimos valores foram assumidos pelo artigo (K2-4)
alfaa = 0.46 // constante
alfac = 0.50 // constante
Kdec = 8.33*10^-4 // valor assumido pelo artigo
q = 1/4*2.5*10^-6  // vazão volumétrica - m³/h
V = 5.5*10^-5 // m³ - volume
F = 96485 // Constante de faraday - C/mol
A = 5*10^-4 // área - m²
YxA = 0.05 // rendimento 
Ecell = 0.70 // volts
Rcell = 4.9389*10^-60 // m³/S /////////////////////////////////////////// n seria m2 em vez de m3 aqui??
M=0.05 //quantidade de microorganismos no início (segundo artigo,0) mol/m³
espce=0.000023 // assumido pelo artigo - m
D =0.08*10^-4*3600
ilim=0.40

solucao=list();
ACET = []
O2 =   []
VOLT = []
CORR = []
ETAA=  []
ETAC = []
POT =  []

original=[q,Rcell,kk,T]
params = ['Vazão q','Rcell','fator preexp kk','Temperatura']
condi = [0,q/2,q*2,Rcell/100,Rcell*100,kk/2,kk*2,T-20,T+20] //Condições para a analise de sensibilidade
for aux=1:size(condi,'c')
    if aux==2 ||aux==3 ; q=condi(aux); Rcell=original(2); kk=original(3);T=original(4) //so alfa muda
    elseif aux==4 ||aux==5 ; Rcell=condi(aux); q=original(1); kk=original(3);T=original(4) //so R muda
    elseif aux==6 || aux==7; kk=condi(aux); q=original(1); Rcell=original(2);T=original(4)// so kk muda
    elseif aux==8 || aux==9 T=condi(aux);q=original(1); Rcell=original(2); kk=original(3) //so T muda
    else q=original(1); Rcell=original(2); kk=original(3); T=original(4) //Nada muda
    end
    
    voltagens = []
    correntes=  []
    etaA = []
    etaC = []
    for vcell=0.454:-0.05:0.18
    
        // Parâmetros do método numérico:
        tolfv = 1.0d-6;  //tolerância na função [mol/tempo]
        tolv = 1.0d-6;  // tolerância na variação de v [mol/volume]
        Itmax = 1000; // número máximo de iterações
        
        k = 0; // número da iteração
        h = 1.0d-8; // incremento usado no cálculo da derivada numérica
        
        ok = 0; // variável que indica critério de parada do while
        
        //getf('funcx.sci'); //declarando função 
        //getd
        
        Jac = zeros(N,N); //dimensionando a matriz Jacobiana 
        
        //Newton-Raphson
        while ok == 0 
          k=k+1;
          fv = funcv(v,v0,kk,K,alfaa,alfac,q,V,F,A,YxA,fx,espce,Kdec,Rcell,vcell,Ecell,M,D,T);    //calcula o vetor fv:
        
          if sum(abs(fv)) <= tolfv then   // checando convergência em f(v)
            ok = 1;
            disp('fv <=tolfv');
          else
            // Calcula a matriz Jacobiana
            for j=1:16
              vj = v(j); // guardando o valor de v(j)
              v(j) = v(j)+h; //incrementando o v(j) em h 
              fvh = funcv(v,v0,kk,K,alfaa,alfac,q,V,F,A,YxA,fx,espce,Kdec,Rcell,vcell,Ecell,M,D,T);
              Jac(:,j)=(fvh-fv)/h;
              v(j) = vj;
            end
            g = Jac\(-fv); // resolução de sistema linear
            v = g + v; //x[k+1]
            
            if sum(abs(g)) <= tolv then  // checando convergência em v
              ok = 2;
              disp('g = delta(v) <= tolv');
            elseif k >= Itmax then  // checando número de iterações
              ok = 3;
              disp('k>=Itmax');
            end
          end
        end
       solucao($+1) = v;
       //salva p curva de desempenho
       etaA = [etaA, v(14)]
       etaC = [etaC, v(15)]
       voltagens = [voltagens,vcell]
       correntes = [correntes,v(16)]
    end
    
    ACET = [ACET,v(1:8)]
    O2 =   [O2, v(9:12)]
    VOLT = [VOLT,voltagens']
    CORR = [CORR,correntes']
    ETAA=  [ETAA, etaA']
    ETAC = [ETAC, etaC']
    POT =  [VOLT.*CORR,(voltagens.*correntes)']
end
// imprindo solução aproximada obtida
//disp('solução')
//v
//disp(v)
// imprimindo o valor da função calculada na solução obtida
//disp('checando função')
//fv = funcv(v,v0,kk,K,alfaa,alfac,q,V,F,A,YxA,fx,espce,Kdec, Rcell, vcell, Ecell, M);
//disp(fv)

disp(solucao) //A lista solucao armazena em cada índice seu um vetor com todas as variáveis de cada iteração de vcell; dessa forma,
              //é possível acessar os valores posteriormente para plotar gráficos, se lembrando da notação do vetor v definido anteriormente.
    
//// Gráficos

alfaa=original(1)
Rcell=original(2)

for aux=0:size(original,'c')-1 //n de parametros q variam na sensibilidade (alfa e RCell)
    index = [1,2+2*aux,3+2*aux]
    t = linspace(1,8,8)
    scf(1+6*aux); plot(t,ACET(t,index))
    legenda = [string(original(1+aux)),string(condi(2+2*aux)),string(condi(3+2*aux))]
    legend(legenda)
    xlabel('y_i')
    ylabel('[Acetato] mol/L')
    
    t2 = linspace(1,4,4) + 8
//    scf(2+6*aux); plot(t2,O2(:,index))
//    legenda = [string(original(1+aux)),string(condi(2+2*aux)),string(condi(3+2*aux))]
//    legend(legenda)
//    xlabel('y_i')
//    ylabel('[O2] mol/L')
    
//    scf(3+6*aux)
    scf(3); subplot(2,2,aux+1);plot(CORR(:,index),VOLT(:,index),Iexp,Vexp,'o')
    legenda = [string(original(1+aux)),string(condi(2+2*aux)),string(condi(3+2*aux)),'exp']
    title(params(aux+1))
    legend(legenda)
    xlabel('Icell A/m2' )
    ylabel('Vcell V')

//    scf(3); subplot(2,2,aux+1);plot(CORR(:,index),VOLT(:,index))
//    legenda = [string(original(1+aux)),string(condi(2+2*aux)),string(condi(3+2*aux))]
//    title(params(aux+1))
//    legend(legenda)
//    xlabel('Icell A/m2' )
//    ylabel('Vcell V')
    
    
//    scf(4); subplot(2,2,aux+1);plot(CORR(:,index),POT(:,index))
//    legenda = [string(original(1+aux)),string(condi(2+2*aux)),string(condi(3+2*aux))]
//    title(params(aux+1))
//    legend(legenda,-4)
//    xlabel('Icell A/m2' )
//    ylabel('Pcell W/m2')
//    
//    scf(5+6*aux); plot(CORR(:,index),ETAA(:,index),[0,12],[0,0])
//    legenda = [string(original(1+aux)),string(condi(2+2*aux)),string(condi(3+2*aux))]
//    legend(legenda)
//    xlabel('Icell A/m2' )
//    ylabel('Sobrepotencial V')
//    title('Anódico')
//
//    scf(6+6*aux); plot(CORR(:,index),ETAC(:,index),[0,12],[0,0])
//    legenda = [string(original(1+aux)),string(condi(2+2*aux)),string(condi(3+2*aux))]
//    legend(legenda)
//    xlabel('Icell A/m2' )
//    ylabel('Sobrepotencial V')
//    title('Catódico')
end

// exporta dados
csvWrite(ACET,'acetato.csv',' ',',')
csvWrite(O2,'oxigenio.csv',' ',',')
csvWrite(VOLT,'voltagens.csv',' ',',')
csvWrite(CORR,'correntes.csv',' ',',')
csvWrite(POT,'potencias.csv',' ',',')
csvWrite(ETAA,'sobrepotencialA.csv',' ',',')
csvWrite(ETAC,'sobrepotencialC.csv',' ',',')







