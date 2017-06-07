% Reproduce Figure 4 from the following paper:
%  
%   Hermes, Nguyen and Winawer (2017). Neuronal synchrony and the relation
%   between the BOLD signal and the local field potential. PLOS Biology
%   http://dx.doi.org/...
%
% Run:
% ns_script07D_Fig4
%
% Purpose: load simulated neural data - time varying membrane potentials
% and show how timeseries is created

%% get all the correlation values 
clear all
sim_nr      = 2;
prm_set_nr  = 1;
elec_nr     = 1;
load(fullfile(boldlfp_RootPath, 'data', ...
    sprintf('NS_simnr%d_elec%d_NS_prmset%d', ...
    sim_nr, elec_nr, prm_set_nr)),'NS');

%% 
numNeurons2plot     = 5;
numNeuronsColor     = {};
trial2plot          = 1;
t_int               = 250:500;

figure('Position',[0 0 600 400]);
ax = get(gca); % to get colororder to make bars the same color as the neuronal timeseries
color_plot = gray(numNeurons2plot+1);
color_plot = color_plot(1:end-1,:);
color_plot = color_plot(end:-1:1,:);

subplot(3,6,1),hold on
for k=1:numNeurons2plot
    y = NS.data.bb_inputs(t_int,k,trial2plot);
    plot(k+zeros(size(y)),'k');
    plot(k+y,'Color',color_plot(k,:))
end
xlim([1 length(t_int)]),ylim([0 6.5])

subplot(3,6,2),hold on
for k=1:numNeurons2plot
    y = NS.data.g_inputs(t_int,k,trial2plot);
    plot(k+zeros(size(y)),'k');
    plot(k+10*y,'Color',color_plot(k,:))
end
xlim([1 length(t_int)]),ylim([0 6.5])

subplot(3,6,3),hold on
for k=1:numNeurons2plot
    y = NS.data.a_inputs(t_int,k,trial2plot);
    plot(k+zeros(size(y)),'k');
    plot(k+8*y,'Color',color_plot(k,:))
end
xlim([1 length(t_int)]),ylim([0 6.5])

subplot(3,6,4),hold on
signal2plot = NS.data.bb_inputs+NS.data.g_inputs+NS.data.a_inputs;
for k=1:numNeurons2plot
    y = signal2plot(t_int,k,trial2plot);
    plot(k+zeros(size(y)),'k');
    plot(k+y,'Color',color_plot(k,:))
end
xlim([1 length(t_int)]),ylim([0 6.5])

subplot(3,6,5),hold on
signal2plot = NS.data.ts(:,:,NS.params.trials_save_inputs(trial2plot));
for k=1:numNeurons2plot
    y = signal2plot(t_int,k);
    plot(k+zeros(size(y)),'k');
    plot(k+y,'Color',color_plot(k,:))
end
xlim([1 length(t_int)]),ylim([0 6.5])

subplot(3,6,6),hold on
signal2plot = NS.data.ts(:,:,NS.params.trials_save_inputs(trial2plot));
for k=1:numNeurons2plot
    barh(k,mean(signal2plot(:,k).^2,1),'FaceColor',color_plot(k,:))
end
ylim([0 numNeurons2plot+1])

subplot(3,6,11),hold on
signal2plot = NS.data.ts(:,:,NS.params.trials_save_inputs(trial2plot));
plot(sum(signal2plot(t_int,:),2),'k')
axis tight

subplot(3,6,12),hold on
signal2plot = NS.data.ts(:,:,NS.params.trials_save_inputs(trial2plot));
barh(1,sum(mean(signal2plot.^2,1),2),'k')
ylim([0 numNeurons2plot+1])

set(gcf,'PaperPositionMode','auto')
print('-dpng','-r300',fullfile(boldlfp_RootPath, 'figures', 'Fig4'))
print('-depsc','-r300',fullfile(boldlfp_RootPath, 'figures', 'Fig4'))


