clear all; close all; clc;

%% ---------Range 45kHz - 60kHz---------%%
wp1=[45e3 60e3]*2*pi;
ws1=[44.5e3 60.5e3]*2*pi;
[N1,wn1]= ellipord(wp1, ws1, 1, 45,'s');
[num1,den1]=ellip(N1, 1, 45, wn1,'s');
sys = tf(num1,den1);

%Decomposing in second-order sections
[sos,g] = tf2sos(num1,den1);
[linhas,colunas] = size(sos);

numx = zeros(linhas,colunas/2);
denx = zeros(linhas,colunas/2);
Sx=g;

for i=1:linhas
    for j=1:colunas/2
        numx(i,j) = sos(i,j);
        denx(i,j) = sos(i,j+3);
    end
    Sx = Sx*tf(numx(i,:),denx(i,:));
end

% In this example we divided the gain equally. However, to implement the real filters it may 
% be better to have more gain in the initial stages. As g<1, we can attenuate the signal 
% putting the gain in the initial stages, avoiding saturation in stages with peaks of gain.

ganho = g^(1/8);

tf1 = ganho*tf(numx(1,:),denx(1,:));
tf2 = ganho*tf(numx(2,:),denx(2,:));
tf3 = ganho*tf(numx(3,:),denx(3,:));
tf4 = ganho*tf(numx(4,:),denx(4,:));
tf5 = ganho*tf(numx(5,:),denx(5,:));
tf6 = ganho*tf(numx(6,:),denx(6,:));
tf7 = ganho*tf(numx(7,:),denx(7,:));
tf8 = ganho*tf(numx(8,:),denx(8,:));

Sx
% figure (1); bodemag(sys,'k'); hold on; bodemag(Sx,'--r'); % Plotting and comparing.

%% -------- Finding the transconductance values -------- %%

% CHOOSING THE SECOND ORDER SECTION 
% It is better to do it manually for each section, as we can control the results changing the arbitrary values. 
% Usually, the values we control are 'g1' and 'g2', but we can also change the capacitances, 'gar' and 'gfr'.

num = tf8.Numerator{1,1};
den = tf8.Denominator{1,1};


%%% -- ARBITRARY VALUES -- %%%

% Reference transconductances 
gar = 1e-5;
gfr = 1e-5;

% Filter order  
N = 2; 
% We use second-order sections, but this topology also allows us to use higher-order sections. 
% If you want to try 4th-order sections, use N = 4 and transform the 8 transfer functions into 4:
% new_tf1 = tf1 * tf2; 
% new_tf2 = tf3 * tf4;
% new_tf3 = tf5 * tf6;
% new_tf4 = tf7 * tf8;
% Then, add the two additional integrators (g3, C3, g4, C4) in the steps below


% Integrators transconductances 
g1 = 1e-7; g2 = 1e-6;  

gx = [g1 g2];

% Integrators capacitances 
C1 = 1e-12; C2 = 1e-12;  

Cx = [C1 C2];

% Integrators - time constants
ctx = [C1/g1 C2/g2];


%%% -- USING THE EQUATIONS FOR THE FLF GM-C FILTER -- %%%

% Calculating Bn
Bn = prod(ctx);

% Numerator Coeficients
Ai = zeros(1, N+1);
for j=1:N+1
        Ai(j) = Bn * num(N+2-j);
end

% Denominator Coeficients
Bi = zeros(1, N+1);
for jj=1:N+1
        Bi(jj) = Bn * den(N+2-jj);
end

% Upper Transconductances
fi = zeros(1, N);
gfi = zeros(1, N);
vet_prod = [];
for k=1:N-1
        vet_prod = ctx(k+1:N);
        fi(k) = Bi(N+1-k)/prod(vet_prod);
        gfi(k) = fi(k)*gfr;
end
fi(N) = Bi(1);
gfi(N) = fi(N)*gfr;

% Lower Transconductances
ai = zeros(1, N+1);
gai = zeros(1, N+1);
ai(1) = Ai(N+1)/Bi(N+1);        % a0 = An/Bn;
ai(N+1) = Ai(1) - ai(1)*fi(N);  % an = A0 - a0*fn;
for kk=1:N-1
        vet_ = ctx(N+1-kk:N);
        ai(N+1-kk) = (Ai(kk+1) - ai(1)*fi(N-kk)*prod(vet_prod))/prod(vet_prod);
end
for kk=1:N+1
        gai(kk) = ai(kk)*gar;
end
