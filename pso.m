%% Inicialize Parpool
clear;
poolobj = gcp;

if isempty(poolobj)
    parpool(parcluster(parallel.defaultClusterProfile));
end
while (isempty(poolobj))
    poolobj = gcp;
end
poolsize = poolobj.NumWorkers;

%% Inicialize parameters
K = 16;           % No. da população
N  = 300;   % No. de passos

cg = 1.5;          % PSO parâmetro CG 
cp = 1.5;          % PSO parâmetro CP

exchange_factor = 0.25;

%% Inicialize Variables
Rp = rand(8,K, N+1);
Rg = rand(8,K, N+1);

f=NaN(K,N+1);
phi=NaN(K,N+1);

x=NaN(8, K, N+1);
x_mean=NaN(8,N+1);

sucess=NaN(K,1);

%Desempenho
ITAE=NaN(2,K,N+1);
Mp=NaN(8,K,N+1);
tss=NaN(4,K,N+1);
d=NaN(N+1,1);

%Display
ITAE_save=NaN(2,N+1);
Mp_save=NaN(8,N+1);
tss_save=NaN(4,N+1);

phi_save=NaN(N+1,1);
phi_pbest_mean=NaN(N+1,1);
phi_mean=NaN(N+1,1);


f_save=NaN(N+1,1);
f_pbest_mean=NaN(N+1,1);
f_mean=NaN(N+1,1);
f_factible = NaN(K,N+1);

g_save=NaN(N+1,1);
x_save=NaN(8,N+1);

vel_abs=NaN(N+1,1);

w_save=NaN(N+1,1);

phi_zero=NaN(N+1,1);

next_exchange_save=NaN(N+1,1);
exchange_save=NaN(N+1,1);


%% Good Solutions
x_good(:,1)=[33.1141;7663.67;9808.72;87.2551;2.0396e-02;-1.2264e-01;5.5544e-01;2.7792e-01];%f=[2072.74;290.429]Mp=[54.68,51.96\39.07,70.76|114.5,47.7\46.4,2.316]tss=[1.71,1.45,8.13,8.22]ITAE=[1932.98,2072.74]
x_good(:,2)=[1.01103;929.598;6841.64;194.23;1.2326e-02;-8.7006e-03;5.6817e-01;3.7716e-01];%f=[1975.22;282.799]Mp=[49.86,48.53\30.99,70.74|109.1,35.97\68.25,8.787]tss=[2.1,1.45,7.84,8.71]ITAE=[1879.76,1975.22]
x_good(:,3)=[1.00959;896.908;8451.91;204.314;1.2485e-02;-9.1335e-03;6.7346e-01;3.7478e-01];%f=[1904.33;274.345]Mp=[49.92,52.19\31.77,70.75|104,43.06\52.57,9.266]tss=[1.75,1.43,8,8.79]ITAE=[1860.25,1904.33]
x_good(:,4)=[1.01113;914.952;6901.21;193.261;1.2333e-02;-8.5724e-03;5.6812e-01;3.6276e-01];%f=[1932.99;273.731]Mp=[53.12,41.56\31.12,69.26|102.7,45.26\61.51,8.573]tss=[1.74,1.28,7.88,8.71]ITAE=[1901.28,1932.99]
x_good(:,5)=[1.01581;934.234;6864.74;181.711;1.2286e-02;-8.9505e-03;5.6112e-01;3.5815e-01];%f=[1930.19;253.874]Mp=[43.3,46.61\31.26,69.91|102.4,29.56\61.76,8.455]tss=[1.59,1.38,7.88,8.71]ITAE=[1901.86,1930.19]
x_good(:,6)=[1.35941;3014.49;9591.72;1.87121;2.5308e-02;-3.8740e-02;6.9031e-01;4.0928e-01];%f=[2031.43;250.125]Mp=[59.13,32.65\33.81,70.75|79.97,45.7\58.93,6.347]tss=[1.83,1.18,7.96,8.52]ITAE=[1963.17,2031.43]
x_good(:,7)=[1.00978;923.278;7954.76;198.228;1.2685e-02;-9.1439e-03;6.4393e-01;3.7691e-01];%f=[1893.71;209.568]Mp=[27.75,34.41\31.45,70.75|86.86,22.32\55.85,8.983]tss=[1.7,1.18,7.96,8.72]ITAE=[1838.25,1893.71]
x_good(:,8)=[6111.11;2.39814;-0.205245;7680.69;1.0858e-01;-1.9348e-01;4.3259e-01;5.5945e-01];%f=[1914.66;129.065]Mp=[56.51,27.45\44.21,45.34|18.56,29.39\12.14,25.37]tss=[1.47,1.15,1.94,9.1]ITAE=[1914.66,1582.23]
x_good(:,9)=[852.506;-39.7808;17.656;150.944;2.3572e-03;1.2451e-04;-1.1394e-04;2.1851e-02];%f=[2172.59;127.474]Mp=[69.56,23.88\59.19,39.37|26.76,26.96\2.627,0.9243]tss=[1.33,1.56,8.38,8.26]ITAE=[2172.59,1614.09]
x_good(:,10)=[229.594;23.9886;-67.6325;159.627;5.9113e-02;-1.3091e-04;1.7724e-04;2.8705e-01];%f=[1442.14;0]Mp=[27.83,25.35\4.113,3.981|17.9,14.52\3.672,2.161]tss=[1.95,2.27,1.32,2.35]ITAE=[1426.61,1442.14]
%x_good(:,11)=[90.6132;5.33609;-14.5871;20.308;2.4295e-02;3.8603e-04;-5.5664e-05;3.7664e-02];%f=[1437.95;0]Mp=[28.87,23.22\4.083,3.979|17.66,15.07\1.984,3.643]tss=[1.76,1.98,1.35,2.37]ITAE=[1437.95,1422.16]
x_good(:,11)=[313.657;23.4986;-90.3365;176.582;0.0797766;-0.000204111;0.000121292;0.27008];%f=[1.4375e+03;0] Mp=[27.93,20.47\4.11,3.94|16.13,15.77\2.502,2.889] ITAE=[1.437e+03;1.412e+03] tss=[2.22;2.34;1.32;2.21]

close all;
figureFullScreen('Name','Results','NumberTitle','off','Color','white');
drawnow;
%% PSO - Main loop
for k=1:N+1
    fprintf('I=%d',k-1);
    %% Update velocity
    v_overflow = 0;
    if k==1
        %% Initialize Position
        x(1,:,1) = abs(normrnd(1000,250,1,K));
        x(2:3,:,1) = normrnd(0,25,2,K);
        x(4,:,1) = abs(normrnd(1000,250,1,K));
        
        x(5,:,1) = abs(normrnd(1e-2,1e-2,1,K));
        x(6:7,:,1) = normrnd(0,1e-4,2,K);
        x(8,:,1) = abs(normrnd(1e-2,1e-2,1,K));
        
        %% Initialize velocity
        vel = zeros(8,K);
        vel(1,:) = 125*rand(1, K);
        vel(2:3,:) = 12.5*rand(2, K);
        vel(4,:) = 125*rand(1, K);
        vel(5,:) = 1e-3*rand(1, K);
        vel(6:7,:) = 5e-5*rand(2, K);
        vel(8,:) = 1e-3*rand(1, K);
        
        %% Initialize exchange index
        exchange=1;
        next_exchange=ceil(exchange_factor*(1+2*rand)*(N-k+2)/(3*(length(x_good)-exchange+1)));
    else
        vel=w*vel+cp*Rp(:,:,k).*(x_pbest-x(:,:,k-1))+cg*Rg(:,:,k).*(repmat(x_gbest,1,K)-x(:,:,k-1));
    end
    %% Update position
    x_overflow = 0;
    if k>1
        x(:,:,k) = x(:,:,k-1) + vel;
        for c=1:8
            if  c==1 || c==4
                upperbound=1e4;
                lowerbound=150;
            elseif c==2 || c==3
                upperbound=100;
                lowerbound=-100;
            elseif c==6 || c==7
                upperbound=1e-3;
                lowerbound=-1e-3;
            else
                upperbound=1;
                lowerbound=0;
            end
            for cparticle=1:K
                loop_count=1;
                while x(c,cparticle,k) > upperbound || x(c,cparticle,k) < lowerbound
                    x_overflow =x_overflow +1;
                    c_j=randi([1 K],2,1);
                    while c_j(1) == c_j(2) || c_j(1) == cparticle || c_j(2) == cparticle
                        c_j=randi([1 K],2,1);
                    end
                    x(c,cparticle,k) = (x_pbest(c,c_j(1))+x_pbest(c,c_j(1))+x_pbest(c,cparticle)+x_gbest(c))/4;
                    if loop_count > 100
                        x(c,cparticle,k) = sign(x(c,cparticle,k))*upperbound*rand;
                        loop_count = 0;
                    end
                    loop_count=loop_count+1;
                end
            end
        end
    end
    
    %% Particle Exchange
    if next_exchange == 0 && exchange <= length(x_good)
        %% Find worst particle
        [worst,worst_I] = max(phi(:,k));
        if worst == 0
            [worst,worst_I] = max(f(:,k));
        end
        x(:,worst_I,k) = x_good(:,exchange);
        fprintf('#');
        %% Update Exchange status
        next_exchange=ceil(exchange_factor*(1+2*rand)*(N-k+2)/(3*(length(x_good)-exchange+1)));
        exchange = exchange + 1;
    end
    
    %% Update particles
    for i = 1:poolsize:K-1,
        parfor j=0:poolsize-1
            [f(i+j,k),phi(i+j,k),ITAE(:,i+j,k),Mp(:,i+j,k),tss(:,i+j,k),sucess(i+j)] = fitness(x(:,i+j,k));
        end
        fprintf('.');
    end

    %% Best particles
    phicount=1;
    if k==1
        %% Pbest
        x_pbest = x(:,:,1);
        f_pbest  = f(:,1)*Inf;
        phi_pbest  = phi(:,1);
        %% Gbest
        if (min(phi_pbest)==0)
            [f_gbest,g] = min(f_pbest);
            phi_gbest = 0;
        else
            f_gbest = min(f_pbest);
            [phi_gbest,g] = min(phi_pbest);
        end
        x_gbest = x_pbest(:,g);
        %% Best
        Mp_best=Mp(:,g,1);
        ITAE_best=ITAE(:,g,1);
        tss_best=tss(:,g,1);
    else
        for i = 1 : K
            %% Pbest
            if phi(i,k) < phi_pbest(i)
                phi_pbest(i) = phi(i,k);
                x_pbest(:,i) = x(:,i,k); 
                phicount=phicount+1;
            end
            if phi(i,k) == 0
                if f(i,k) < f_pbest(i)
                    f_pbest(i)  = f(i,k); 
                    x_pbest(:,i) = x(:,i,k);
                end
            end
            %% ALL - Gbest
            for j = 1:k         
                if phi(i,j) < phi_gbest
                    %% Gbest - Phi
                    phi_gbest = phi(i,j);
                    x_gbest = x(:,i,j);
                    %% Best - Phi
                    Mp_best = Mp(:,i,j);
                    tss_best = tss(:,i,j);
                    ITAE_best = ITAE(:,i,j);
                    %f_gbest = sum(ITAE_best);
                end
                if phi(i,j) == 0
                    if f(i,j) < f_gbest
                        %% Gbest - f
                        f_gbest = f(i,j);
                        x_gbest = x(:,i,j);
                        %% Best - f
                        Mp_best = Mp(:,i,j);
                        tss_best = tss(:,i,j);
                        ITAE_best = ITAE(:,i,j);
                    end
                end
            end
        end
    end
    %% Update Inertia
    d(k)=max(std(x_pbest,0,2,'omitnan'));
    w=0.9-0.4*d(k)/max(d,[],'omitnan');
    
    %% Save results
    ITAE_save(:,k)=ITAE_best;
    Mp_save(:,k)=Mp_best;
    tss_save(:,k)=tss_best;

    phi_save(k)=phi_gbest;
    phi_pbest_mean(k)=mean(phi_pbest);
    phi_mean(k)=mean(phi(:,k));
    
    f_save(k)=f_gbest;
    f_pbest_mean(k)=mean(f_pbest(~isinf(f_pbest)));
    
    g_save(k)=max(ITAE_best);
    
    x_save(:,k)=x_gbest;
    
    vel_abs(k,:) = mean(sqrt(sum(vel.^2,2)));
    
    x_mean(:,k)=mean(x(:,:,k),2,'omitnan');
    
    w_save(k)=0.9-0.4*d(k)/max(d,[],'omitnan');
    next_exchange_save(k)=next_exchange-1;
    exchange_save(k)=exchange;
    
    phi_zero(k)=sum(phi(:,k)==0);
    
    %% Plot results
    
        %% Evolution of parameters
        %ITAE
        subaxis(4,4,1,'MR',0.01,'ML',0.05,'MT',0,'MB',0.05,'ML',0.05);
        semilogy(0:k-1,ITAE_save(1,1:k),'r.-');
        hold on;
        semilogy(0:k-1,ITAE_save(2,1:k),'b.-');
        hold off;
        ylabel('ITAE');
        legend('ITAE^d','ITAE^q','Orientation','horizontal','Location','northoutside');
        grid;
        axis(gca,'tight');

        %% tss
        %tss SPS
        subaxis(4,4,2,'MR',0.01,'ML',0.05,'MT',0,'MB',0.05,'MT',0,'MB',0.05);
        semilogy(0:k-1,tss_save(1,1:k),'r.-');
        hold on;
        semilogy(0:k-1,tss_save(2,1:k),'b.-');
        hold off;
        ylabel('t_{ss} SPS');
        legend('t_{ss}^d','t_{ss}^q','Orientation','horizontal','Location','northoutside');
        grid;
        axis(gca,'tight');
        
        %tss SS
        subaxis(4,4,3,'MR',0.01,'ML',0.05,'MT',0,'MB',0.05,'MT',0,'MB',0.05);
        semilogy(0:k-1,tss_save(3,1:k),'g.-');
        hold on;
        semilogy(0:k-1,tss_save(4,1:k),'m.-');
        hold off;
        ylabel('t_{ss} SS');
        legend('t_{ss}^d','t_{ss}^q','Orientation','horizontal','Location','northoutside');
        grid;
        axis(gca,'tight');
        
        %% Mp
        
        %Mp SPS
        subaxis(4,4,4,'MR',0.01,'ML',0.05,'MT',0,'MB',0.05);
        semilogy(0:k-1,Mp_save(1,1:k),'r.-');
        hold on;
        semilogy(0:k-1,Mp_save(2,1:k),'b.-');
        semilogy(0:k-1,Mp_save(5,1:k),'g.-');
        semilogy(0:k-1,Mp_save(6,1:k),'m.-');
        hold off;

        ylabel('M_p SPS');
        legend('M_{p,d}','M_{p,q}','M_{p,d}^{coup}','M_{p,q}^{coup}','Orientation','horizontal','Location','northoutside');
        grid;
        axis(gca,'tight');
        
        %Mp SS
        subaxis(4,4,5,'MR',0.01,'ML',0.05,'MT',0,'MB',0.05);
        semilogy(0:k-1,Mp_save(3,1:k),'r.-');
        hold on;
        semilogy(0:k-1,Mp_save(4,1:k),'b.-');
        semilogy(0:k-1,Mp_save(7,1:k),'g.-');
        semilogy(0:k-1,Mp_save(8,1:k),'m.-');
        hold off;

        ylabel('M_p SS');
        legend('M_{p,d}','M_{p,q}','M_{p,d}^{coup}','M_{p,q}^{coup}','Orientation','horizontal','Location','northoutside');
        grid;
        axis(gca,'tight');

        %% Fitness
        %f
        subaxis(4,4,6,'MR',0.01,'ML',0.05,'MT',0,'MB',0.05);
        semilogy(0:k-1,g_save(1:k),'b.-');
        hold on;
        semilogy(0:k-1,f_save(1:k),'m.-');
        semilogy(0:k-1,f_pbest_mean(1:k),'g.-');
        hold off;
        axis(gca,'tight');
        
        ylabel('Unbounded Fitness');
        legend('g_{g}','f_{g}','f_{p}','Orientation','horizontal','Location','northoutside');
        grid;

        %phi
        subaxis(4,4,7,'MR',0.01,'ML',0.05,'MT',0,'MB',0.05);
        semilogy(0:k-1,phi_save(1:k),'r.-');
        hold on;
        semilogy(0:k-1,phi_pbest_mean(1:k),'b.-');
        semilogy(0:k-1,phi_mean(1:k),'m.-');
        hold off;
        ylabel('Bounded fitness');
        legend('\phi_g','\phi_p','\phi','Orientation','horizontal','Location','northoutside');
        grid;
        axis(gca,'tight');

        %% Evolution of abs(Vel)
        
        subaxis(4,4,8,'MR',0.01,'ML',0.05,'MT',0,'MB',0.05);
        semilogy(0:k-1,vel_abs(1:k),'r.-');
        ylabel('Mean of |v|_2');
        legend('mean(|v|_2)','Orientation','horizontal','Location','northoutside');
        grid;
        axis(gca,'tight');
        
        %% No. of particles outside of restrictions
        subaxis(4,4,9,'MR',0.01,'ML',0.05,'MT',0,'MB',0.05);
        bar(0:k-1,phi_zero(1:k),'c','EdgeColor','none');
        ylabel('Feasiable particles');
        legend('No. of phi(x)=0','Orientation','horizontal','Location','northoutside');
        grid;
        xlim([0 max(k-1,1)]);
        ylim([0 max(1,min(inf,max(phi_zero(1:k))))]);
        
        %% Evolution of Ex
        subaxis(4,4,10,'MR',0.01,'ML',0.05,'MT',0,'MB',0.05);
        yyaxis left;
        plot(0:k-1,next_exchange_save(1:k),'r.-');
        axis(gca,'tight');
        ylabel('Next');
        set(gca,'ycolor','r');
        yyaxis right;
        plot(0:k-1,exchange_save(1:k)-1,'b.-');
        axis(gca,'tight');
        ylabel('Done');
        set(gca,'ycolor','b');
        legend('Next','Done','Orientation','horizontal','Location','northoutside');
        grid;

        %% Evolution of W
        %Wy
        subaxis(4,4,11,'MR',0.01,'ML',0.05,'MT',0,'MB',0.05);
        plot(0:k-1,x_save(1,1:k),'r.-');
        hold on;
        plot(0:k-1,x_save(2,1:k),'b.-');
        plot(0:k-1,x_save(3,1:k),'g.-');
        plot(0:k-1,x_save(4,1:k),'m.-');
        hold off;
        legend('W_{y,11}','W_{y,12}','W_{y,21}','W_{y,22}','Orientation','horizontal','Location','northoutside');
        grid;
        axis(gca,'tight');
        
        %Wu
        subaxis(4,4,12,'MR',0.01,'ML',0.05,'MT',0,'MB',0.05);
        plot(0:k-1,x_save(5,1:k),'r.-');
        hold on;
        plot(0:k-1,x_save(6,1:k),'b.-');
        plot(0:k-1,x_save(7,1:k),'g.-');
        plot(0:k-1,x_save(8,1:k),'m.-');
        hold off;

        ylabel('W_u');
        legend('W_{u,11}','W_{u,12}','W_{u,21}','W_{u,22}','Orientation','horizontal','Location','northoutside');
        grid;
        axis(gca,'tight');
        
        %% Evolution of W_mean
        subaxis(4,4,13,'MR',0.01,'ML',0.05,'MT',0,'MB',0.05);
        plot(0:k-1,x_mean(1,1:k),'r.-');
        hold on;
        plot(0:k-1,x_mean(2,1:k),'b.-');
        plot(0:k-1,x_mean(3,1:k),'g.-');
        plot(0:k-1,x_mean(4,1:k),'m.-');
        hold off;

        ylabel('Mean of W_y');
        legend('W_{y,11}','W_{y,12}','W_{y,21}','W_{y,22}','Orientation','horizontal','Location','northoutside');
        grid;
        axis(gca,'tight');

        subaxis(4,4,14,'MR',0.01,'ML',0.05,'MT',0,'MB',0.05);
        plot(0:k-1,x_mean(5,1:k),'r.-');
        hold on;
        plot(0:k-1,x_mean(6,1:k),'b.-');
        plot(0:k-1,x_mean(7,1:k),'g.-');
        plot(0:k-1,x_mean(8,1:k),'m.-');
        hold off;

        ylabel('Mean of W_u');
        legend('W_{u,11}','W_{u,12}','W_{u,21}','W_{u,22}','Orientation','horizontal','Location','northoutside');
        grid;
        axis(gca,'tight');

        %% Evolution of d
        subaxis(4,4,15,'MR',0.01,'ML',0.05,'MT',0,'MB',0.05);
        semilogy(0:k-1,d(1:k),'r.-');

        ylabel('Evolution of d');
        legend('d','Orientation','horizontal','Location','northoutside');
        grid;
        axis(gca,'tight');

        %% Evolution of inertia
        subaxis(4,4,16,'MR',0.01,'ML',0.05,'MT',0,'MB',0.05);
        plot(0:k-1,w_save(1:k),'b.-');

        ylabel('Inertia');
        legend('Inertia','Orientation','horizontal','Location','northoutside');
        grid;
        axis(gca,'tight');
        
        drawnow;
    %% Print last Results
    fprintf('/%.4g|%d|%f/x_o=%d|v_o=%d/S=%d|Next=%d|Ex=%d\n',w,phicount,d(k),x_overflow,v_overflow,sum(sucess),next_exchange,exchange-1);
    fprintf('Wy=[%g,%g;%g,%g] Wu=[%.4e,%.4e;%.4e,%.4e]\n',x_gbest);
    fprintf('x=[%g;%g;%g;%g;%g;%g;%g;%g]\n',x_gbest);
    fprintf('f=[%.4e;%g] Mp=[%.4g,%.4g\\%.4g,%.4g|%.4g,%.4g\\%.4g,%.4g] ITAE=[%.3e;%.3e] tss=[%.4g;%.4g;%.4g;%.4g]\n\n',f_gbest,phi_gbest,Mp_best,ITAE_best,tss_best);
    
    %% Update Exchanges
    if next_exchange >0
        if k>2 && phi_save(k-1) == phi_gbest && f_save(k) == f_save(k-1)
            next_exchange = next_exchange - 1;
        elseif exchange <= length(x_good)
            next_exchange=ceil(exchange_factor*(1+2*rand)*(N-k+2)/(3*(length(x_good)-exchange+1)));
        end
    else
        next_exchange=0;
    end
end
fprintf('\n');