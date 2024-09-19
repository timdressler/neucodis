% pp_wPLI_analysis.m
%
% Calculation and plotting of phase coherence as indexed by the wPLI (Vinck et al., 2015) 
%
% Originally by Dr. Mike Cohen
% Adapted from Tim Dressler, 19.09.2024

clear
close all
clc

%select electrodes for the analysis
electrodes2use = { 'Fz';'Pz' };

%convert channel labels to indices
elecs2use = zeros(size(electrodes2use));
for i=1:length(electrodes2use)
    elecs2use(i) = find(strcmpi(electrodes2use{i},{EEG.chanlocs.labels}));
end

%seeded phase synchronization (set to empty ("{}") for no analyses)
electrodes4seeded_synch = { chanlocs(elecs2use(1)).labels , chanlocs(elecs2use(2)).labels };

%wavelet parameters
min_freq =  2;
max_freq = 40;
num_frex = 25;
wavelet_cycle_range = [ 3 12 ];

frex = logspace(log10(min_freq),log10(max_freq),num_frex);


%gaussian width and time
s = logspace(log10(wavelet_cycle_range(1)),log10(wavelet_cycle_range(2)),num_frex)./(2*pi.*frex);
t = -2:1/EEG.srate:2;

%fft and convolution details
Ltapr  =  length(t);
Ldata  =  prod(ntime*ntrials);
Lconv1 =  Ldata+Ltapr-1;
Lconv  =  pow2(nextpow2( Lconv1 ));

wavelets = zeros(num_frex,length(t));
for fi=1:num_frex
    wavelets(fi,:) = exp(2*1i*pi*frex(fi).*t).*exp(-t.^2./(2*s(fi)^2));
end

%initialize output and hidden-layer TF matrices
tf    = zeros(nchans+2,num_frex,length(times2saveidx),2);
allphasevals = zeros(nchans,num_frex,length(times2save),ntrials,2);
synchOverTrials = zeros(2,length(electrodes4seeded_synch),nchans,num_frex,length(times2saveidx),2);
allAS = zeros(2,num_frex,ntime,ntrials,2);

%run convolution
% loop around channels
for chani=1:nchans+2
    
    %fft of data (channel or true dipoles)
    if chani<=nchans
        EEGfft = fft(reshape(EEG.data(chani,:,:),1,[]),Lconv);
        Lapfft = fft(reshape(simulatedLap(chani,:,:),1,[]),Lconv);
    else
        [EEGfft,Lapfft] = deal(fft(reshape(sourceTimeSeries(:,mod(chani,nchans),:),1,[]),Lconv));
    end
        
    
    %loop over frequencies and complete convolution
    for fi=1:num_frex
        
        %average reference
        %convolve and get analytic signal (as)
        as = ifft(EEGfft.*fft(wavelets(fi,:),Lconv),Lconv);
        as = as(1:Lconv1);
        as = reshape(as(floor((Ltapr-1)/2):end-1-ceil((Ltapr-1)/2)),ntime,ntrials);
        
        %enter into TF matrix
        temppow = mean(abs(as).^2,2);
        tf(chani,fi,:,1) = 10*log10( temppow(times2saveidx)/mean(temppow(baseidx(1):baseidx(2))) );
        
        %save phase values
        if chani<=nchans
            allphasevals(chani,fi,:,:,1) = as(times2saveidx,:);
        end
        
        %all values from all time points
        if chani==elecs2use(1)
            allAS(1,fi,:,:,1) = as;
        elseif chani==elecs2use(2)
            allAS(2,fi,:,:,1) = as;
        end
        
        %Laplacian
        %convolve and get analytic signal (as)
        as = ifft(Lapfft.*fft(wavelets(fi,:),Lconv),Lconv);
        as = as(1:Lconv1);
        as = reshape(as(floor((Ltapr-1)/2):end-1-ceil((Ltapr-1)/2)),ntime,ntrials);
        
        %enter into TF matrix
        temppow = mean(abs(as).^2,2);
        tf(chani,fi,:,2) = 10*log10( temppow(times2saveidx)/mean(temppow(baseidx(1):baseidx(2))) );
        
        %save phase values
        if chani<=nchans
            allphasevals(chani,fi,:,:,2) = as(times2saveidx,:);
        end
        
        %all values from all time points
        if chani==elecs2use(1)
            allAS(1,fi,:,:,2) = as;
        elseif chani==elecs2use(2)
            allAS(2,fi,:,:,2) = as;
        end
        
    end %end frequency loop
end %end channel loop

%compute phase connectivity over trials

for chanx=1:length(electrodes4seeded_synch)
    
    %ISPC (average reference and laplacian)
    synchOverTrials(1,chanx,:,:,:,1) = mean(exp(1i* bsxfun(@minus,angle(allphasevals(strcmpi(electrodes4seeded_synch{chanx},{chanlocs.labels}),:,:,:,1)),angle(allphasevals(:,:,:,:,1))) ),4);
    synchOverTrials(1,chanx,:,:,:,2) = mean(exp(1i* bsxfun(@minus,angle(allphasevals(strcmpi(electrodes4seeded_synch{chanx},{chanlocs.labels}),:,:,:,2)),angle(allphasevals(:,:,:,:,2))) ),4);
    
    
    %wPLI (average reference and laplacian)
    cdd = bsxfun(@times,allphasevals(strcmpi(electrodes4seeded_synch{chanx},{chanlocs.labels}),:,:,:,1),conj(allphasevals(:,:,:,:,1)));
    cdi = imag(cdd);
    synchOverTrials(2,chanx,:,:,:,1) = mean( abs( mean( abs(cdi).*sign(cdi) ,4) )./mean(abs(cdi),4) ,4);
    
    cdd = bsxfun(@times,allphasevals(strcmpi(electrodes4seeded_synch{chanx},{chanlocs.labels}),:,:,:,2),conj(allphasevals(:,:,:,:,2)));
    cdi = imag(cdd);
    synchOverTrials(2,chanx,:,:,:,2) = mean( abs( mean( abs(cdi).*sign(cdi) ,4) )./mean(abs(cdi),4) ,4);
end

%NaN's cause electrode data shifts in the eeglab topoplot function
synchOverTrials(isnan(synchOverTrials))=0;

%% plot data

times2plot = dsearchn(times2save',[200 800]');
freq2plot  = dsearchn(frex',centfreq);

clim = .8;

figure(1), clf, colormap hot
subplot(221)
topoplot(squeeze(mean(mean(abs(synchOverTrials(1,1,:,freq2plot-1:freq2plot+1,times2plot(1):times2plot(2),1)),5),4)),chanlocs,'maplimits',[0 clim],'plotrad',.63,'numcontour',0,'style','map','electrodes','off','emarker2',{elecs2use '.' 'g' 8 1});
title('seeded connectivity, average reference')

subplot(222)
topoplot(squeeze(mean(mean(abs(synchOverTrials(1,1,:,freq2plot-1:freq2plot+1,times2plot(1):times2plot(2),2)),5),4)),chanlocs,'maplimits',[0 clim],'plotrad',.63,'numcontour',0,'style','map','electrodes','off','emarker2',{elecs2use '.' 'g' 8 1});
title('seeded connectivity, Laplacian')

subplot(223)
contourf(times2save,frex,abs(squeeze(synchOverTrials(1,1,elecs2use(2),:,:,1))),40,'linecolor','none')
set(gca,'clim',[0 clim])
hold on
toplot = squeeze(mean(abs(synchOverTrials(1,1,elecs2use(2),freq2plot-1:freq2plot+1,:,1)),4));
plot(times2save,toplot*30+10,'w','linew',2)
plot(get(gca,'xlim'),[10 10],'w--')

subplot(224)
contourf(times2save,frex,abs(squeeze(synchOverTrials(1,1,elecs2use(2),:,:,2))),40,'linecolor','none')
set(gca,'clim',[0 clim])
hold on
toplot = squeeze(mean(abs(synchOverTrials(1,1,elecs2use(2),freq2plot-1:freq2plot+1,:,2)),4));
plot(times2save,toplot*30+10,'w','linew',2)
plot(get(gca,'xlim'),[10 10],'w--')
xlabel('Time (ms)'), ylabel('Frequency (Hz)')


figure(2), clf, colormap hot
subplot(221)
topoplot(squeeze(mean(mean(abs(synchOverTrials(2,1,:,freq2plot-1:freq2plot+1,times2plot(1):times2plot(2),1)),5),4)),chanlocs,'maplimits',[0 clim],'plotrad',.63,'numcontour',0,'style','map','electrodes','off','emarker2',{elecs2use '.' 'g' 8 1});
title('seeded connectivity, average reference')

subplot(222)
topoplot(squeeze(mean(mean(abs(synchOverTrials(2,1,:,freq2plot-1:freq2plot+1,times2plot(1):times2plot(2),2)),5),4)),chanlocs,'maplimits',[0 clim],'plotrad',.63,'numcontour',0,'style','map','electrodes','off','emarker2',{elecs2use '.' 'g' 8 1});
title('seeded connectivity, Laplacian')

subplot(223)
contourf(times2save,frex,abs(squeeze(synchOverTrials(2,1,elecs2use(2),:,:,1))),40,'linecolor','none')
hold on
toplot = squeeze(mean(abs(synchOverTrials(2,1,elecs2use(2),freq2plot-1:freq2plot+1,:,1)),4));
plot(times2save,toplot*30+10,'w','linew',2)
plot(get(gca,'xlim'),[10 10],'w--')
set(gca,'clim',[0 clim])

subplot(224)
contourf(times2save,frex,abs(squeeze(synchOverTrials(2,1,elecs2use(2),:,:,2))),40,'linecolor','none')
hold on
toplot = squeeze(mean(abs(synchOverTrials(2,1,elecs2use(2),freq2plot-1:freq2plot+1,:,1)),4));
plot(times2save,toplot*30+10,'w','linew',2)
plot(get(gca,'xlim'),[10 10],'w--')
set(gca,'clim',[0 clim])
xlabel('Time (ms)'), ylabel('Frequency (Hz)')

 