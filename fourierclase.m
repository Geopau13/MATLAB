clear all
close all
clear clf
clc
dy = 1/12;
yi = 1953;
ye = 1963;
X = yi:dy:ye-dy;
X = X';
data = importdata('recife.dat');

YO = data;
N = length(YO);
M = N/2; %se determinan la cantidad dea aronicos totales a obtener
t = [1:N]'; % vectro del tiempo o vectro de candiad de datos
f0 = 1/N; % frecuencia fundamental
fN = M/N; % frecuendcia de Nuquist

for k = 0:M;
    ak(k+1) = (2/N)*sum(YO.*cos(t*2*pi*f0*k)); %se obtienen las amplitudes AK
    bk(k+1) = (2/N)*sum(YO.*sin(t*2*pi*f0*k)); %se obtienen las amplitudes BK
end


    
    ak(1) = ak(1)/2; %este es el primer amronico a0    

    % No nos interesa el armónico fundamental
a0 = ak(1); %se guarda a0
ak = ak(2:end); %se reescribe el vector que almacena los ak
bk = bk(2:end); %se reescribe el vector que almacena los bk

% Calculamos las amplitudes y las fases
ck = sqrt (ak.*ak + bk.*bk); % pitagoras
fk = atan2(bk,ak)*180/pi; % se tiene las fases ya en grados
freq = [(1:N/2)/N]'; %se obtiene el vector de frec, hasta N/2.

figure('Units', 'Normalized',...
'Position',[0.1 0.13 0.5 0.4],...
'PaperPositionMode','auto',...
'color','w')
plot(X,YO,'linewidth',2)
set(gca,'xlim',[1953 1963],'ylim',[22 30],'ytick',[22:1:30],...
'ticklength',[0.02 04],'fontsize',15,...
'yminortick','on','xminortick','on','xtick',[1953:1:1963])
xlabel('Tiempo (Años)','fontsize',16)
ylabel('Tempratura (°C)','fontsize',16)

figure('Units', 'Normalized',...
'Position',[0.1 0.13 0.5 0.4],...
'PaperPositionMode','auto',...
'color','w')
plot(freq,ck,'-o','linewidth',1.5,...
'markeredgecolor','b','markerfacecolor','r')
set(gca,'xlim',[0 0.5],'ylim',[0 1.6],'ytick',[0:0.5:1.6],...
'ticklength',[0.02 04],'fontsize',15,...
'yminortick','on','xminortick','on','xtick',[0:0.1:0.5])
xlabel('Frecuencia (cph)','fontsize',15)
ylabel('Amplitud (°C)','fontsize',15)

T = 1./freq; % se obtienen los periodos
varex = 100*ck.^2/sum(ck.^2); % procentaje de la varizna explicada
matriz=[[1:M]' T freq ak' bk' ck' fk' varex'];

fprintf('---------------------------------------------------------------\n')
fprintf('Armo Periodo Frec ak bk ck fk Var\n')
fprintf(' k mes cpm °C °C °C Grados Percent\n')
fprintf('---------------------------------------------------------------\n')
fprintf(' %2d %7.3f %5.3f %6.3f %6.3f %7.3f %8.3f %6.3f\n',matriz')

