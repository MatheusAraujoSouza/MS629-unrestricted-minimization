function x=projeto629_21112020
clc;
%Matheus araujo Souza ra 184145 , ms629 projeto 1

%Primeira parte do algoritmo vamos construir nossas
%funções de direção de descida  
%%%%%%%%%%%%%%%%%%ponto inicial%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
x=[-22];
%%%%%%%%%%%%%%%%% declando quais sao as variaveis simbolicas%%%%%%%%%%%%%%%
syms x1



%%%%%%%%%%%%%%%%%%%espaço para colocar as funcoes %%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%Quadratica%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%fun =1*x1^2+ 2*x2^2+3*x3^2;
%fun =1*x1^2 + 2*x2^2;
%%%%%%%%%%%%%%%%%%%%%%%%Styblinsky-tang%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%fun = x1^4 - 16*x1^2 + 5*x1 
%fun = x1^4 - 16*x1^2 + 5*x1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%Rosenbrook%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%fun=100*(x2-x1^2)^2 +(x1 - 1)^2;
%%%%%%%%%%%%%%%%%%%%%%%%%%%Rastring%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fun=x1^2-10*cos(2*pi*x1);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Quando estiver aumentando ou diminuindo as variaveis lembrar de mudar os
%x1 ou x2 para x3 .... ou apenas x1 respectivamente, logo lembrar sempre de
%ir dendro das funcoes e mudar os valores dos xi quando for mudar o tamanho
%do n 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
gra=gradient(fun,[x1]);
He=hessian(fun,[x1]);

%%%%%%%%%%%%%%%%%%%espaço para invocar as funções%%%%%%%%%%%%%%%%%%%%%%%%%%
%tic
%Gradiente(fun,gra,x);
%toc
tic
Newton(fun,gra,He,x)
toc
tic
DFP(fun,gra,x)
toc
tic
Gradiente(fun,gra,x);
toc
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%sobre comandos usados%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%alguns comando usados que devem ser levados em consideracao 
%1)subs e' usado para trazer de volta as variaveis syms para valores type
%float de double precison
%2)vpa (vpa(x)usa aritmética de ponto flutuante de precisão variável (VPA) 
%para avaliar cada elemento da entrada simbólica xem pelo menos ddígitos 
%significativos, onde dé o valor da digitsfunção.O valor padrão digitsé 32.
%3)norm(v) retorna a norma euclidiana do vetor v.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%funcao%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%Gradiente%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function G=Gradiente(fun,gra,x)
syms x1 


alpha=1.0*10^-4;
gamma=0.5;
epsilon=10^-6;
k=0;

M=1000; % numero maximo de iteracoes.

while ((vpa(norm(subs(gra,[x1],vpa(x))),13) >= epsilon) & (k < M) )
     
       
        d= - vpa(subs(gra,[x1],vpa(x)),13);
        t=1; %passo inicial
        while(vpa(subs(fun,[x1],vpa(x + t*d)),13) > vpa(subs(fun,[x1],vpa(x)),13) + vpa(alpha*t*(subs(gra,[x1],vpa(x))')*d,13))
          t=t*gamma;
        end
        x = x + t * d;
        k=k+1;
        
    
end

Vk=k
Vx=vpa(x,5)
VF=vpa(subs(fun,[x1],vpa(x,5)),5)
VG=vpa(subs(gra,[x1],vpa(x,5)),5)

    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%Newton%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function N=Newton(fun,gra,He,x)
syms x1 
sigma=10^-3;
alpha=1.0*10^-4;
gamma=0.5;
teta=10^-6;
Beta= 10^-3;
epsilon=10^-6;
k=0;

[m,n]=size(subs(gra,[x1],vpa(x)));
I=eye(m);


M=1000; % numero maximo de iteracoes.

while ((norm(subs(gra,[x1],vpa(x))) >= epsilon) & (k < M) )
    mi=0;
    R=vpa(chol(subs(He,[x1],vpa(x)) + mi*I),4);
    d=R\(R'\(-subs(gra,[x1],vpa(x))));
    while(subs(gra,[x1],vpa(x))'* d >= - teta* norm(subs(gra,[x1],vpa(x)))*norm(d))
        mi=max(2*mi,Beta);
        R=vpa(chol(subs(He,[x1],vpa(x)) + mi*I),4);
        d=R\(R'\(-subs(gra,[x1],vpa(x))));
    end
    if(norm(d) < sigma*norm(subs(gra,[x1],vpa(x))))
        d = ((sigma*norm(subs(gra,[x1],vpa(x))))/norm(d))*d;
    end
    t=1;
    while(subs(fun,[x1],vpa(x + t*d)) > subs(fun,[x1],vpa(x)) + alpha*t*subs(gra,[x1],vpa(x))'*d)
          t= t*gamma;
    end
    x = x + t * d;
    k=k+1;
        
end
Vk=k
Vx=vpa(x,5)
VF=vpa(subs(fun,[x1],vpa(x,5)),5)
VG=vpa(subs(gra,[x1],vpa(x,5)),5)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%DFP%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function D=DFP(fun,gra,x)
%%%%%%%%%Aqui não vamos precisar usar a Hessiana%%%%%%%%%%%%%%%%%%%%%%%%%%%
syms x1 

sigma=10^-3;
alpha=10^-4;
gamma=0.5;
teta=10^-6;
epsilon=10^-6;
Beta= 10^-3;


%agora vamos definir o valor de Ho
[m,n]=size(subs(gra,[x1],vpa(x)));
Ho=eye(m);
k=0;
M=1000; % numero maximo de iteracoes.
while ((vpa((norm(subs(gra,[x1],vpa(x)))),13) >= epsilon) & (k < M) )
    d = - vpa(Ho*subs(gra,[x1],vpa(x)),13);
    if((subs(gra,[x1],vpa(x))')* d >= - teta* norm(subs(gra,[x1],vpa(x)))*norm(d))
        d= - vpa(subs(gra,[x1],vpa(x)),13);
        Ho=eye(m);
    end
    
    if(norm(d) < vpa(sigma*norm(subs(gra,[x1],vpa(x))),13))
        d = vpa(((sigma*norm(subs(gra,[x1],vpa(x))))/norm(d))*d,13);
    end
    
    t=1;
    while(vpa(subs(fun,[x1],vpa(x + t*d)),13) > vpa(subs(fun,[x1],vpa(x)),13) + vpa(alpha*t*(subs(gra,[x1],vpa(x))')*d,13))
          t= gamma*t;
    end
    
    tempx=x;
     x = x + t*d;
     p = x - tempx;
     q = vpa(subs(gra,[x1],vpa(x)),13) - vpa(subs(gra,[x1],vpa(tempx)),13);
     r=Ho*q;
     
     if((p')*q > 0 )
         Ho = Ho + (p*(p'))/((p') * q) - (r * (r'))/((r') * q);
     end
   k=k+1;
  
end
Vk=k
Vx=vpa(x,5)
VF=vpa(subs(fun,[x1],vpa(x,5)),5)
VG=vpa(subs(gra,[x1],vpa(x,5)),5)















