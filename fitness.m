function [f,phi,ITAE,Mp,tss,sucess,Signal] = fitness(x,benchmark,delta)
% Simulations parameters
switch nargin
    case 0
        benchmark=1;
        x=[30.8589;589.029;1834.55;474.089;9.0483e-04;3.3655e-02;-4.6392e-02;3.5266e-04];
        x=[1e2;0;0;1e2;1e-02;0;0;1e-02];
        x=[852.506;-39.7808;17.656;150.944;0.00235724;0.000124509;-0.000113942;0.0218507];
        x=[850;-30;30;150;0.002;0;0;0.02];
        x=[1;0;0;1;0.05/400;0;0;0.6/400];
        x=[1;0;0;1;1e-3;0;0;1e-3];
        x=[1;0;0;1;1e-8;0;0;1e-8];
        x=[365.157;-299.161;212.965;742.962;9.0922e-03;-2.5575e-02;-1.6894e-02;9.5126e-03];
        delta=[0,0];
        x=[1;0;0;1;1e-3;0;0;1e-3];
        x=[1;0;0;1;1;0;0;1];
        x=[1;0;0;1;10^0.45;0;0;10^0.45];
        x=[313.657;23.4986;-90.3365;176.582;0.0797766;-0.000204111;0.000121292;0.27008];
    case 1
        delta=[0,0];
        benchmark=0;
    case 2
        delta=[0,0];
    otherwise
end
MpSS=10;
MpSSq=5;
MpSPS=35;
MpSPSq=25;

Ny=20; %Horizonte de predição
Nu=10; %Horizonte de controle
Wyp = zeros(2*Ny,2*Ny);
Wup = zeros(2*Nu,2*Nu);
T=1e-4;
Variac=120; %Tensão no barramento

dLm=delta(1);
dRr=delta(2);

%% Step Values

StepInit=[1,1]; %Valor inicial da resposta ao degrau
StepFinal=[3,3]; %Valor final da resposta ao degrau

%StepInit=[5.10,2.25];
%StepFinal=[4.32,1.13];

StepTime=[1,1.5]; %Tempo em que ocorre a resposta ao degrau
timeaxis=[0.75,1.75];
%% Weighting matrices
Wy=[x(1),x(2);x(3),x(4)];
Wu=[x(5),x(6);x(7),x(8)];

for i=1:Ny
    Wyp((i-1)*2+1:i*2,(i-1)*2+1:i*2) = Wy;
end

for i=1:Nu
    Wup((i-1)*2+1:i*2,(i-1)*2+1:i*2) = Wu;
end

%% Simulation
try
    simOut = sim('dfig','SrcWorkspace','Current');
    time = simOut.get('tout');
    Erro = [simOut.get('Erro'),simOut.get('ErroIdeal')];
    Sinal = [simOut.get('Sinal'),simOut.get('SinalIdeal')];
    SizeWindow = [0.004,0.49,0.49];
    sucess=1; %Simulação ocorreu com sucesso
catch
    sucess=0;
end

%% Find parameters
tss=64*[1;1;1;1]; %tempo de acomodação
ITAE=100000*[1;1]; %Indice de desempenho
Mp=10000*[1;1;1;1;1;1;1;1]; %Sobressinal

if sucess==1
    if benchmark == 1
        close all;
        figure('Name','Time Response','NumberTitle','off');
        drawnow;
        subaxis(1,2,1,'MR',0.01,'ML',0.05,'MT',0.02,'MB',0.05); %
        I_cuttime=[find(time>timeaxis(1),1);find(time>timeaxis(2),1)];
        Signal.time=time(I_cuttime(1):8:I_cuttime(2)); 
        Signal.SPSd=Sinal(I_cuttime(1):8:I_cuttime(2),1);
        Signal.SPSq=Sinal(I_cuttime(1):8:I_cuttime(2),2);
        Signal.SSd=Sinal(I_cuttime(1):8:I_cuttime(2),3);
        Signal.SSq=Sinal(I_cuttime(1):8:I_cuttime(2),4);
    end
    %% RealModel
    I_StepTime=[0,0];
    NoiseMax=[0,0];
    NoiseStd=[0,0];
    Mp_I=[0,0];
    Mp_v=[0,0];
    tss_I=[0,0];
    for WindowIncrease = 0:0.06:0.06
        I_TimeWindow=zeros(2);
        for i=1:2
            I_StepTime(i)=find(time>StepTime(i),1);
            I_TimeWindow(i,1)=find(time>StepTime(i)+SizeWindow(1)+WindowIncrease,1,'first');
            I_TimeWindow(i,2)=find(time>StepTime(i)+SizeWindow(2),1,'first');
            
            temp=abs(Erro(I_TimeWindow(i,1):I_TimeWindow(i,2),i));
            NoiseMax(i)=max(temp);
            NoiseStd(i)=std(temp);
            temp = Sinal(I_StepTime(i):I_TimeWindow(i,1),i)-StepFinal(i);
            temp = abs(temp);
            temp_2 = abs(Erro(I_StepTime(i):I_TimeWindow(i,2),i));
            ITAE(i)=sum((time(I_StepTime(i):I_TimeWindow(i,2))-time(I_StepTime(i))).*temp_2);%max(temp_2,temp_2.^2)
            
            [temp_2,Mp_I(i)] = max(Sinal(I_StepTime(i):I_TimeWindow(i,1),i));
            Mp(i) = abs(abs(temp_2-StepFinal(i))/(StepInit(i)-StepFinal(i)))*100;
            Mp_v(i)=temp_2;
            temp = I_StepTime(i)+find(temp>=NoiseMax(i)+NoiseStd(i),1,'last');
            temp=temp+find(Sinal(temp:I_TimeWindow(i,1),i)-StepFinal(i)<0,1);
            temp=temp+find(Sinal(temp:I_TimeWindow(i,1),i)-StepFinal(i)>0,1);
            
            if ~isempty(temp)
                tss_I(i)=temp;
                tss(i)=(time(temp)-StepTime(i))*1000;
            end
        end
        temp=abs(Erro(I_StepTime(2):I_TimeWindow(2,2),1));
        temp=sum((time(I_StepTime(2):I_TimeWindow(2,2))-time(I_StepTime(2))).*temp);%max(temp,temp.^2)
        ITAE(1) = ITAE(1)+temp;
        temp=abs(Erro(I_StepTime(1):I_TimeWindow(1,2),2));
        temp=sum((time(I_StepTime(1):I_TimeWindow(1,2))-time(I_StepTime(1))).*temp);%max(temp,temp.^2)
        ITAE(2) = ITAE(2)+temp;
        temp_2 = max(abs(Erro(I_StepTime(2):I_TimeWindow(2,2),1)));
        Mp(5)=abs(temp_2/(StepInit(2)-StepFinal(2))*100);
        temp_2 = max(abs(Erro(I_StepTime(1):I_TimeWindow(1,2),2)));
        Mp(6)=abs(temp_2/(StepInit(1)-StepFinal(1))*100);
        if abs(tss_I)~=0
            break;
        end
    end
    if benchmark == 1
        for i=1:2
            if  tss_I(i)==0
                tss_I(i)=I_TimeWindow(i,1);
            end
            subplot(1,2,i); %
            temp=[linspace(StepInit(i),StepInit(i),length(time(1:I_StepTime(i)))),linspace(StepFinal(i),StepFinal(i),length(time(I_StepTime(i)+1:end)))]; %
            plot(time,Sinal(:,i),'r') %
            hold on; %
            plot(time,temp,'k') %
            scatter(time(I_StepTime(i)+Mp_I(i)),Mp_v(i),'g','d'); %
            scatter(time(tss_I(i)),Sinal(tss_I(i),i),'c'); %
            plotWindow=[StepTime(i)-2/1000,StepTime(i)+WindowIncrease+0.00401];
            xlim(plotWindow);
            if i==1
                ylabel('Current d axis (A)');
                if tss(i)==64
                    title(sprintf('t_{ss,d}>=%g ms|M_{p,d}=%.4g %%',tss(i),Mp(i)));
                else
                    title(sprintf('t_{ss,d}=%g ms|M_{p,d}=%.4g %%',tss(i),Mp(i)));
                end
            else
                if tss(i)==64
                    title(sprintf('t_{ss,q}>=% gms|M_{p,q}=%.4g %%',tss(i),Mp(i)));
                else
                    title(sprintf('t_{ss,q}=% gms|M_{p,q}=%.4g %%',tss(i),Mp(i)));
                end
                ylabel('Current q axis (A)');
            end
            legend('Signal','Reference','Estimated overshoot','Estimated settling time','Location','best')
            xlabel('Time(s)');
            grid;
        end
    end
    %% IdealModel
    NoiseMax=[0,0];
    NoiseStd=[0,0];
%     if benchmark == 1
%         subaxis(1,2,2,'MR',0.01,'ML',0.05,'MT',0.02,'MB',0.05); %
%     end
    for i=1:2
        I_TimeWindow(i,1)=find(time>StepTime(i)+SizeWindow(1)+0.06,1,'first');
        I_TimeWindow(i,2)=find(time>StepTime(i)+SizeWindow(2),1,'first');
        temp=abs(Erro(I_TimeWindow(i,1):I_TimeWindow(i,2),i+2));
        NoiseMax(i)=max(temp);
        NoiseStd(i)=std(temp);
        temp_2 = max(0.02*(StepFinal(i)-StepInit(i)),NoiseMax(i)+NoiseStd(i)); %NoiseMax(i)+NoiseStd(i);
%         if benchmark == 1
%             temp=[linspace(StepInit(i),StepInit(i),length(time(1:I_StepTime(i)))),linspace(StepFinal(i),StepFinal(i),length(time(I_StepTime(i)+1:end)))]; %
%             plot(time,linspace(temp_2,temp_2,length(time))+StepFinal(i),time,Sinal(:,i+2)) %
%             hold on; %
%             plot(time,temp) %
%         end
        temp = Sinal(I_StepTime(i):I_TimeWindow(i,2),i+2)-StepFinal(i);
        
        temp = I_StepTime(i)+find(temp>=temp_2,1,'last');
        temp_3 = find(Sinal(temp:I_TimeWindow(i,1),i+2)-StepFinal(i)<0,1);
        if ~isempty(temp_3)
            temp=temp+temp_3;
        end
        temp_3 = find(Sinal(temp:I_TimeWindow(i,1),i+2)-StepFinal(i)<0,1);
        if ~isempty(temp_3)
            temp=temp+temp_3;
        end
        
        [temp_2,temp_I] = max(Sinal(I_StepTime(i):I_TimeWindow(i,2),i+2));
        Mp(i+2) = abs(abs(temp_2-StepFinal(i))/(StepInit(i)-StepFinal(i)))*100;
        
        if ~isempty(temp)
            tss(i+2)=(time(temp)-StepTime(i))*1000;
%             if benchmark == 1
%                 scatter(time(I_StepTime(i)+temp_I),temp_2); %
%                 scatter(time(temp),Sinal(temp,i+2)); %
%             end
        end
    end
    temp_2 = max(abs(Erro(I_StepTime(2):I_TimeWindow(2,2),3)));
    Mp(7)=abs(temp_2/(StepInit(2)-StepFinal(2))*100);
    temp_2 = max(abs(Erro(I_StepTime(1):I_TimeWindow(1,2),4)));
    Mp(8)=abs(temp_2/(StepInit(1)-StepFinal(1))*100);
end

%%Result
f = max(ITAE); %função objetivo
phi = 0; %função de restrições
for i=1:2
    if Mp(i) > MpSPS
        phi = phi + Mp(i)-MpSPS;
    end
    if Mp(i+2) > MpSS
        phi = phi + Mp(i+2)-MpSS;
    end
    if Mp(i+4) > MpSPSq
        phi = phi + Mp(i+4)-MpSPSq;
    end
    if Mp(i+6) > MpSSq
        phi = phi + Mp(i+6)-MpSSq;
    end
    temp = max(tss(i),tss(i+2));
    if temp > 3
        phi = phi + temp-3;
    end
end
if benchmark==1
    fprintf('Mp=[%.4g,%.4g\\%.4g,%.4g|%.4g,%.4g\\%.4g,%.4g] tss=[%g,%g,%g,%g]\n',Mp,tss);
    fprintf('f=[%.4e,%.4g] ITAE=[%.4e,%.4e] Wy=[%g,%g;%g,%g] Wu=[%g,%g;%g,%g]\n',f,phi,ITAE,Wy,Wu);
end