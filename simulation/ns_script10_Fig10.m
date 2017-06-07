%%
% This script runs the simulation used in Hermes et al:
%
%
% DH 2016

clear all

%% load data and simulation and correlate across all frequencies

% load ECoG/fMRI data structure:
load(fullfile(boldlfp_RootPath, 'data', 'boldecog_structure_final'));

v_area = zeros(length(data),1);

% correlate data across all frequencies
r_data = NaN(length(data),200);
for l = 1:length(data)
    v_area(l) = data{l}.v_area;
    ecog_spectralchange = squeeze(median(data{l}.allboots_ecog250,2));
    bold_change = data{l}.betas;
    r_data(l,:) = corr(ecog_spectralchange,bold_change');
end

% correlate simulation across all frequencies 
r_simulation = NaN(length(data),200);
for l = 1:length(data)
    % load output from the first model (BB - level, G - coh, A - level)
    prm_set = 1;
    load(fullfile(boldlfp_RootPath, 'data', ['NS_simnr' int2str(sim_nr) '_elec' int2str(l) '_NS_prmset' int2str(prm_set)]))

    simlfp_spectralchange = squeeze(median(NS.data.lfp_spectra_bs,2));
    simbold_change = median(NS.data.bold_bs,2);
    r_simulation(l,:) = corr(simlfp_spectralchange,simbold_change);  
end

%% 

figure('Position',[0 0 800 500])  

% ---- DATA: Plot BOLD correlation with ECoG across frequencies -----
subplot(2,2,1), set(gca, 'FontSize', 10),hold on
plot(f,zeros(size(f)),'k','LineWidth',1)
plot(f,r_data(ismember(v_area,1),:),'Color',[.5 .5 .5],'LineWidth',1)
plot(f,r_data(21,:),'r','LineWidth',2)
plot(f,mean(r_data(ismember(v_area,1),:),1),'k','LineWidth',2)
xlabel('Frequency (Hz)'), ylabel('correlation (r)')
xlim([0 200])
ylim([-1 1])
title('A Data V1')

subplot(2,2,2), set(gca, 'FontSize', 10),hold on
plot(f,zeros(size(f)),'k','LineWidth',1)
plot(f,r_data(ismember(v_area,[2 3]),:),'Color',[.5 .5 .5],'LineWidth',1)
plot(f,r_data(18,:),'r','LineWidth',2)
plot(f,mean(r_data(ismember(v_area,[2 3]),:),1),'k','LineWidth',2)
xlabel('Frequency (Hz)'), ylabel('correlation (r)')
xlim([0 200])
ylim([-1 1])
title('B Data V2/V3')


% ---- SIMULATION: Plot BOLD correlation with LFP across frequencies -----
subplot(2,2,3), set(gca, 'FontSize', 10),hold on
plot(f,zeros(size(f)),'k','LineWidth',1)
plot(f,r_simulation(ismember(v_area,1),:),'Color',[.5 .5 .5],'LineWidth',1)
plot(f,r_simulation(21,:),'r','LineWidth',2)
plot(f,mean(r_simulation(ismember(v_area,1),:),1),'k','LineWidth',2)
xlabel('Frequency (Hz)'), ylabel('correlation (r)')
xlim([0 200])
ylim([-1 1])
title('C Simulations fit to V1')

subplot(2,2,4), set(gca, 'FontSize', 10),hold on
plot(f,zeros(size(f)),'k','LineWidth',1)
plot(f,r_simulation(ismember(v_area,[2 3]),:),'Color',[.5 .5 .5],'LineWidth',1)
plot(f,r_simulation(18,:),'r','LineWidth',2)
plot(f,mean(r_simulation(ismember(v_area,[2 3]),:),1),'k','LineWidth',2)
xlabel('Frequency (Hz)'), ylabel('correlation (r)')
xlim([0 200])
ylim([-1 1])
title('D Simulations fit to V2/V3')


set(gcf,'PaperPositionMode','auto')
fname = fullfile(boldlfp_RootPath, 'figures', 'Fig10');
print('-dpng','-r300',fname)
print('-depsc','-r300',fname)
