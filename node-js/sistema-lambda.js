//Desenvolvido por Marco Polo Moreno de Souza em 24/02/2022
//Resolve equações de Bloch no domínio do tempo
//PIBIC Patryky - Passagem adiabática por Raman estimulado: aplicação no relógio atômico de rubídio

fs = require('fs');

const nome_arquivo = '/dados.txt'

const path = __dirname + `${nome_arquivo}`

function lerArquivo(caminho){
    fs.readFile(caminho, 'utf-8', function(error,data){
        if(error){
            console.log('erro de leitura: ' + error.message)
        } else {
            console.log(data)
        }

    })
}

function escreverArquivo(caminho,texto){
    fs.writeFile(caminho, texto, function(error){
        if (error){
            console.error('erro de escrita' + error.message)
        } else {
            console.log('escreve com sucesso em '+ caminho)
        }
    })
}

//lerArquivo(path)



Pi = 3.141592653589793

q33 = 2*Pi*5e6;
q13 = 0.5*q33;
q12 = 0;

omega21 = 2*Pi*6835e6;
w1 = -0.5*omega21;
w2 =  0.5*omega21;

LB = 2*Pi*0.5e12      //Largura de banda, em Hz
LM = 2*Pi*1000e6      //Largura da mordida espectral, em Hz
A = 1.4e11;

d = 0;

h = 50e-15;
pontos = 5e6

t = -h*pontos/2;
a11 = 0.5, a22 = 0.5;
a33 = 0, a12 = 0, b12 = 0;
a13 = 0, b13 = 0, a23 = 0, b23 = 0;

var alpha


function bloch(a11, a22, a33, a12, b12, a13, b13, a23, b23, j)  //sistema de 3 níveis lambda
{
    if (j===1) return 2*Omega*(Math.cos(alpha)*b13 - Math.sin(alpha)*a13)                +0.5*q33*a33;
    if (j===2) return 2*Omega*(Math.cos(alpha)*b23 - Math.sin(alpha)*a23)                +0.5*q33*a33;
    if (j===3) return 2*Omega*(Math.sin(alpha)*(a13+a23) - Math.cos(alpha)*(b13+b23))        -q33*a33;
    
    if (j===4) return -(w2-w1)*b12 + Omega*(Math.cos(alpha)*(b13+b23) - Math.sin(alpha)*(a13+a23)) - a12*q12;
    if (j===5) return  (w2-w1)*a12 + Omega*(Math.cos(alpha)*(a23-a13) + Math.sin(alpha)*(b23-b13)) - b12*q12;
   
    if (j===6) return  -(d-w1)*b13 + Omega*Math.cos(alpha)*b12 + Omega*Math.sin(alpha)*(a12+a11-a33) - a13*q13;
    if (j===7) return   (d-w1)*a13 + Omega*Math.sin(alpha)*b12 + Omega*Math.cos(alpha)*(a33-a11-a12) - b13*q13;
    if (j===8) return  -(d-w2)*b23 - Omega*Math.cos(alpha)*b12 + Omega*Math.sin(alpha)*(a12+a22-a33) - a23*q13;
    if (j===9) return   (d-w2)*a23 - Omega*Math.sin(alpha)*b12 + Omega*Math.cos(alpha)*(a33-a22-a12) - b23*q13;
}

function campo(T)
{
    alpha = -0.5*omega21*T

    return 2*Pi*2*A * ( Math.sin(T*LB) - Math.sin(T*LM) ) / (T*LB);
}

k1 = [], k2 = [], k3 = [], k4 = [];

dados = "tempo rho11 rho22 rho33 sigma12 sigma13 sigma23 soma\n" + "ps\n"

for (k=0; k<=pontos; k++){

    Omega = campo(t)

    for (p=1; p<=9; p++)
        k1[p] = bloch( a11, a22, a33, a12, b12, a13, b13, a23, b23, p );
    
    Omega = campo(t + 0.5*h)

    for (p=1; p<=9; p++)
        k2[p] = bloch( a11 + 0.5*h*k1[1], a22 + 0.5*h*k1[2], a33 + 0.5*h*k1[3], 
                        a12 + 0.5*h*k1[4], b12 + 0.5*h*k1[5], a13 + 0.5*h*k1[6], 
                        b13 + 0.5*h*k1[7], a23 + 0.5*h*k1[8], b23 + 0.5*h*k1[9], p );

    Omega = campo(t + 0.5*h)

    for (p=1; p<=9; p++)
        k3[p] = bloch( a11 + 0.5*h*k2[1], a22 + 0.5*h*k2[2], a33 + 0.5*h*k2[3], 
                        a12 + 0.5*h*k2[4], b12 + 0.5*h*k2[5], a13 + 0.5*h*k2[6], 
                        b13 + 0.5*h*k2[7], a23 + 0.5*h*k2[8], b23 + 0.5*h*k2[9], p );

    Omega = campo(t + h)

    for (p=1; p<=9; p++)
        k4[p] = bloch( a11 + h*k3[1], a22 + h*k3[2], a33 + h*k3[3], 
                        a12 + h*k3[4], b12 + h*k3[5], a13 + h*k3[6], 
                        b13 + h*k3[7], a23 + h*k3[8], b23 + h*k3[9], p );

    t = t + h;
    
    a11 = a11 + h/6.0 * ( k1[1] + 2*k2[1] + 2*k3[1] + k4[1] );
    a22 = a22 + h/6.0 * ( k1[2] + 2*k2[2] + 2*k3[2] + k4[2] );
    a33 = a33 + h/6.0 * ( k1[3] + 2*k2[3] + 2*k3[3] + k4[3] );

    a12 = a12 + h/6.0 * ( k1[4] + 2*k2[4] + 2*k3[4] + k4[4] );
    b12 = b12 + h/6.0 * ( k1[5] + 2*k2[5] + 2*k3[5] + k4[5] );

    a13 = a13 + h/6.0 * ( k1[6] + 2*k2[6] + 2*k3[6] + k4[6] );
    b13 = b13 + h/6.0 * ( k1[7] + 2*k2[7] + 2*k3[7] + k4[7] );

    a23 = a23 + h/6.0 * ( k1[8] + 2*k2[8] + 2*k3[8] + k4[8] );
    b23 = b23 + h/6.0 * ( k1[9] + 2*k2[9] + 2*k3[9] + k4[9] );

    soma = a11 + a22 + a33
    sigma12 = a12*a12 + b12*b12
    sigma13 = a13*a13 + b13*b13
    sigma23 = a23*a23 + b23*b23

    if (k%1000 === 0){
        console.log((1e12*t).toFixed(8) + " " + a11.toFixed(4) + " " + a22.toFixed(4) + " " + a33.toFixed(4) + " " + soma.toFixed(4));
        dados = dados + 1e12*t + " " + a11 + " " + a22 + " " + a33 + " " + sigma12 + " " + sigma13 + " " + sigma23 + " " + soma + " " + campo(t) + "\n"
    }    
}

escreverArquivo(path, dados)