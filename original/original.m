% Created by Eugene M. Izhikevich, February 25, 2003
% modified by E.Lowet, 5th September 2014
Ne=4000;                 Ni=1000;          % number of excitatory and inhibitory neurons
a=[0.02*ones(Ne,1);     0.1*ones(Ni,1)];   % a =tiemscale of recover varibale u
b=[0.2*ones(Ne,1);      0.2*ones(Ni,1)];   % b= sensitivity of  u to subthreshold oscillations
c=[-65*ones(Ne,1);      -65*ones(Ni,1)];   % c= membrane voltage after spike (reset)
d=[8*ones(Ne,1);          2*ones(Ni,1)];   % d= spike reset of recover varibale u
v=-65*ones(Ne+Ni,1);                       % Initial values of v = voltage
u=b.*v;                                    % Initial values of u= membrane recovery variable
firings=[];                                % spike timings
simulation_time=8000 ;
dt=1;Ntot=(Ne+Ni);
disp('defining neurons done')
%%%%%%%%%%%%%%  gaussian input%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
var_E= 1.5;%% to excitatory neurons  
var_I= 1.5;%% to inhibitory neurons
%% creating  main input stimulus
clear stim_input
Amplitude=1; % sinusoidal spatial modualtion of input strength
Meanlevel=7; % mean input level to RS cells
stim_input(1:Ne,1) = (sin((-(1*pi):(2*pi)/((Ne./1)-1):(1*pi)))).*Amplitude+Meanlevel;
stim_input(Ne+1:Ntot,1)= ones(Ni,1).*3.5; % additional mean inputto FS cells
%%%%%%%%%%%%%%%%%%%% synaptic constants %%%%%%%%%%%%%%%%%%%%%
gampa=zeros(Ne,1,'single');
gaba= zeros(Ni,1,'single');
decay_ampa =1;decay_gaba =4;
rise_ampa =0.1;rise_gaba =0.1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Constructing connectivity matrix %%%%%%%%%%%%%%%%%%%%%%%%
disp('start constructing connectivity matrix')
EE = 0.004;see=0.4; %% ecitatory  to excitatory
EI = 0.07;sei=0.3;  %% ecitatory  to inhibitory
IE =-0.04;sie=0.3;%% inhibitory to excitatory
II =-0.015;sii=0.3; %% inhibitory to inhibitory
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
S=zeros(Ntot,'single');dist=zeros(Ne,'single');
op=-pi:(2*pi)/(Ne-1):pi;
for ind=1:Ne
    dist(ind,:) =   abs( angle(   exp(1i*op(ind)  )./exp(1i*op))  );
end
dist((dist<0.001))= NaN; %
KEE= (EE*exp(-dist./see));
%%%%%%%%%%%%%%%%%%%%%%%%%%
dist=zeros(Ni,'single');
op=-pi:(2*pi)/(Ni-1):pi;
for ind=1:Ni
    dist(ind,:) = abs( angle(   exp(1i*op(ind)  )./exp(1i*op))  );
end
dist((dist<0.001))= NaN;
KII= (II*exp(-dist./sii));
%%%%%%%%%%%%%%%%%%%%%%%%%%
dist=zeros(Ne,Ni,'single');
op2=-pi:(2*pi)/(Ne-1):pi;op=-pi:(2*pi)/(Ni-1):pi;
for ind=1:Ne
    dist(ind,:) =   abs( angle(   exp(1i*op2(ind)  )./exp(1i*op))  );
end
dist((dist<0.001))= NaN;
KEI= (EI*exp(-dist./sei));
%%%%%%%%%%%%%%%%%%%%%%
dist=zeros(Ni,Ne,'single');
op2=-pi:(2*pi)/(Ne-1):pi;op=-pi:(2*pi)/(Ni-1):pi;
for ind=1:Ni
    dist(ind,:) =abs( angle(   exp(1i*op(ind)  )./exp(1i*op2))  );
end
dist((dist<0.001))= NaN;
KIE= (IE*exp(-dist./sie));
%%%%%%%%%%%%%%%%%%%%%%
S(1:Ne,1:Ne)= KEE;
S(Ne+1:Ntot,Ne+1:Ntot)= KII;
S(1:Ne,Ne+1:Ntot)= KIE';
S(Ne+1:Ntot,1:Ne)= KEI';
S(isnan(S))=0;S=single(S);
clear dist  KII KIE KEI
disp('done constructing connectivity matrix')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%% !!!!!! MAIN LOOP !!!!!!!!  %%%%%%%%%%%%%%%%%%%%%%%%
disp('start simulation')
for t=1:dt:simulation_time
    if mod(t,25) ==0
    disp([ num2str(t) 'ms of ' num2str(simulation_time ) 'ms'])
    end
     I=[var_E*randn(Ne,1);var_I*randn(Ni,1)]+stim_input; % thalamic input
    fired=find(v>=30);    % indices of spikes
    firings=[firings; t+0*fired,fired];
    v(fired)=c(fired);
    u(fired)=u(fired)+d(fired);
    %synaptic potentials
    gampa=gampa + dt*(0.3*(((1+tanh((v(1:Ne)/10)  +2 ))/2).*(1-gampa)/rise_ampa - gampa/decay_ampa));
    gaba=  gaba + dt*(0.3*(((1+tanh((v(Ne+1:end)/10)  +2 ))/2).*(1-gaba)/rise_gaba - gaba/decay_gaba));
    gsyn=[gampa ;gaba];
    % defining input to eah neuron as the summation of all synaptic input
    % form all connected neurons
    I=I+S*gsyn;
    v=v+0.5*(0.04*v.^2+5*v+140-u+I); % step 0.5 ms
      v=v+0.5*(0.04*v.^2+5*v+140-u+I); % for numerical
    u=u+a.*(b.*v-u);                 % stability
end;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('done simulation')
 
%%%%%%%%%%%% *AFTER SIMULATION ANALYSIS* %%%%%%%%%%%%%%%%%%%%%%%
firings=int32(firings); 
%%%%%%%%%%%%******** Plottting%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('start plotting')
 
figure('Color','w','Position' ,[ 100 100  600 350]),
subplot(2,1,1,'Fontsize',15) % spike  raster
firingexc=firings(find(firings(:,2) <=Ne),:);firinginh=firings(find(firings(:,2) > Ne),:);
plot(firingexc(:,1),firingexc(:,2),'.','Color', [ 0.8 0.2 0.2]); % spike raster
hold on, plot(firinginh(:,1),firinginh(:,2),'.','Color', [ 0.2 0.2 0.8]);
axis tight;set(gca,'xticklabel',[])
xlim([ 600 1800])
subplot(2,1,2,'Fontsize',15) 
Fs = 1000./dt;[t1,t2] = hist(firings(:,1),0:1:t);
spectrogram(((t1)-mean(t1)),252,250,20:0.5:50,Fs,'Yaxis');
axis xy
xlim([ 0.6 1.8]);xlabel('Time s')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% Making a spike matrix and computing spike rate %%%%%%%%%%%%%%%
clear spikerate spikerate2
rastersp=zeros(Ne,max(firingexc(:,1)),'int8');
rastersp2=zeros(Ni,max(firinginh(:,1)),'int8');
nn=0;
for ind=  (1):(Ne) % RS
    nn=nn+1;
    rastersp(nn,firingexc(find(firingexc(:,2)==ind),1).*(1/dt))=1;  
end
nn=0;
for ind=  (Ne+1):(Ne+Ni) % FS
    nn=nn+1;
    rastersp2(nn,firinginh(find(firinginh(:,2)==ind),1).*(1/dt))=1;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('computation cross-correlatio matrix')
%%%%% Creating full phase-locking and phase difference matrix %%%%%%%%
    spike_dat=rastersp2;
    clear allcoh alltim
    nn1=0;timwin=200:simulation_time-50;
    for seed=(1:10:size(rastersp2,1))%steps of 10 to speed up
        disp([num2str(seed) ' of ' num2str(size(rastersp2,1))])
        nn1=nn1+1;nn2=0;
        for ind=(1:10:size(rastersp2,1))%1
            nn2=nn2+1;
            if seed ~= ind
                sig1= (double(spike_dat(seed:seed,timwin))'); %
                sig2= (double(spike_dat(ind:ind,timwin))');
                [c,lags]=xcorr( sig1 ,  sig2    ,12,'coeff'); % here +/- 12ms
                [num1 num2]= max(c);
                allcoh(nn1,nn2)= num1; %peak height
                alltim(nn1,nn2)=num2;  %peak lag
            else  % for the autocorrelation case
                allcoh(nn1,nn2)= 1;
                alltim(nn1,nn2)=13;
            end
        end
    end
    figure('COlor','w','Position',[300 300 240 200]),subplot(1,1,1,'Fontsize',15);
    imagesc(allcoh); % spike cross-correlation peak
    colormap('hot');%colorbar
    set(gca,'CLim', [0 0.4])
       set(gca,'xticklabel',[],'yticklabel',[]);
       figure('COlor','w','Position',[300 300 240 200]),
    phs=((alltim-13).*(-1));
    for ind=1:size(phs,1)
        for ind2=1:size(phs,2)
          phs(ind,ind2)=  phs(ind,ind2)./(  (spikerate2(ind)+spikerate2(ind2))/2   ./2).*pi;%
        end
    end 
    subplot(1,1,1,'Fontsize',17);h=imagesc(phs); % spike timing difference
    acoh=allcoh;
    tt=(acoh)>0.1 &  (acoh)<1;  % phase-locking threshold (here arbitrarly defined)
    set(h,'AlphaData',tt );
    set(h, 'AlphaDataMapping', 'scaled');
    set(gca,'xticklabel',[],'yticklabel',[]);
   set(gca,'Clim',[-pi./2 pi./2])
