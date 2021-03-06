function ns_script03_Fig3
% Reproduce Figure 3 for the following paper:
%  
%   Hermes, Nguyen and Winawer (2017). Neuronal synchrony and the relation
%   between the BOLD signal and the local field potential. PLOS Biology
%   http://dx.doi.org/...
%
% Run:
% ns_script03_Fig3
%
% This script loads the data, plots a power spectrum from the ECoG data and
% shows the broadband, gamma and alpha signal changes during each condition
% for one electrode at the time
%
%
% DH 2017

%%
clear all
close all

load(fullfile(boldlfp_RootPath, 'data','boldecog_structure_final'),'data');

%% plot power spectrum for one electrode

elec = 7; % the electrode to plot

figure('Position',[0 0 150 150]),hold on
for k = 1:length(data{elec}.ecog_spectra250_500) % [1 4 7]
    plot(data{elec}.ecog_f,mean(data{elec}.ecog_spectra250_500(k).data,1),...
        'Color',data{elec}.colors{k},'LineWidth',2)
end

set(gca, 'XScale', 'log', 'YScale', 'log','XTick',[10 100])
set(gca,'YScale', 'log','XTick',[10 100])

xlabel ('Frequency (Hz)'), ylabel('Power')
xlim([5 200])

set(gcf,'PaperPositionMode','auto')
print('-painters','-r300','-dpng',strcat(['./figures/Fig3B_el' int2str(elec)]));
print('-painters','-r300','-depsc',strcat(['./figures/Fig3B_el' int2str(elec)]));

%% plot broadband, gamma and alpha changes

elec = 7;

figure('Position',[0 0 100 150])
subplot(3,1,1),hold on
y = mean(data{elec}.bb_all,2);
y_err = [quantile(data{elec}.bb_all,.16,2) quantile(data{elec}.bb_all,.84,2)];
for k = 1:length(data{elec}.labels)
    bar(k,y(k),'FaceColor',data{elec}.colors{k})
    plot([k k],y_err(k,:),'k')
end
xlim([0 11]),ylim([min([y;0])-.1 max([y;0])+.1])
set(gca,'XTick',[1:1:10],'YTick',[-1:.5:1])

subplot(3,1,2),hold on
y = mean(data{elec}.gamma_all,2);
y_err = [quantile(data{elec}.gamma_all,.16,2) quantile(data{elec}.gamma_all,.84,2)];
for k = 1:length(data{elec}.labels)
    bar(k,y(k),'FaceColor',data{elec}.colors{k})
    plot([k k],y_err(k,:),'k')
end
xlim([0 11]),ylim([min([y;0])-.1 max([y;0])+.1])
set(gca,'XTick',[1:1:10],'YTick',[-1:.5:1])

subplot(3,1,3),hold on
y = mean(data{elec}.alpha_all,2);
y_err = [quantile(data{elec}.alpha_all,.16,2) quantile(data{elec}.alpha_all,.84,2)];
for k = 1:length(data{elec}.labels)
    bar(k,y(k),'FaceColor',data{elec}.colors{k})
    plot([k k],y_err(k,:),'k')
end
xlim([0 11]),ylim([min(y)-.1 max([y;0])+.1])
set(gca,'XTick',[1:1:10],'YTick',[-1:.5:1])

set(gcf,'PaperPositionMode','auto')
print('-painters','-r300','-dpng',strcat(['./figures/Fig3C_el' int2str(elec)]));
print('-painters','-r300','-depsc',strcat(['./figures/Fig3C_el' int2str(elec)]));
