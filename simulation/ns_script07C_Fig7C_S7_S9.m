%%
% This script generates pannels for Fig 7C, Supplementary Figures S7 and S9
% from Hermes et al:
%
% Purpose: load simulated neural data - time varying membrane potentials -
% that are fit to ECoG data and plot results
%
% DH 2016
clear all

%% load simulations fitted to ECoG data from all electrodes

sim_nr = 2;
els = 1:1:22;

%%% OUTPUTS:
v_area = NaN(length(els),1); %visual area per electrode

r2_data_fit   = NaN(8,length(els)); % R2 between BOLD data and fit for each model:
r2_data_fit_s = NaN(8,length(els)); % R2 from scambling predictions
r2_data_data  = NaN(1,length(els));% % R2 from BOLD12 to BOLD34

% DATA: 4:bb,g,a,bold, nr els, up to 10 conditions, 8 simulations
all_data = NaN(4,length(els),10); % ECoG / BOLD data
% SIMULATION: 4:bb,g,a,bold, nr els, up to 10 conditions, 8 simulations
all_simulation = NaN(4,length(els),10,8); % BOLD simulation

% load the ECoG/fMRI data
load(fullfile(boldlfp_RootPath, 'data', 'boldecog_structure_final'));

for l = 1:length(els)
    
    elec = els(l);
       
    % load the simulation outputs 
    load(fullfile(boldlfp_RootPath, 'data', sprintf('NS_simnr%d_elec%d_simulation_outputs', sim_nr, elec)),'simulation_outputs');

    v_area(l) = data{l}.v_area;

    data_bb = median(data{elec}.bb_all,2);
    data_g = median(data{elec}.gamma_all,2);
    data_a = median(data{elec}.alpha_all,2);
    data_bold = data{elec}.betas * mean(data{elec}.norm);
    
    all_data(1,l,1:length(data_bb))=data_bb;
    all_data(2,l,1:length(data_g))=data_g;
    all_data(3,l,1:length(data_a))=data_a;
    all_data(4,l,1:length(data_bold))=data_bold;
    
    % get simulated ECoG (bb, g, a) and BOLD responses into
    % 'all_simulation' matrix
    for k=1:8 
        % ECoG predictions:
        all_simulation(1,l,[1:size(simulation_outputs,1)],k) = simulation_outputs(:,k,1);
        all_simulation(2,l,[1:size(simulation_outputs,1)],k) = simulation_outputs(:,k,2);
        all_simulation(3,l,[1:size(simulation_outputs,1)],k) = simulation_outputs(:,k,3);
        % BOLD, raw from simulation:
        all_simulation(4,l,[1:size(simulation_outputs,1)],k) = simulation_outputs(:,k,4);
        
        fitted_bold = simulation_outputs(:,k,4);        
        r2_data_fit(k,l) = ns_cod(fitted_bold,data_bold', true);
                
        tmp = NaN(1,100);
        for ii = 1:100
            tmp(ii) = ns_cod(fitted_bold(randperm(length(fitted_bold))),data_bold', true);
            r2_data_fit_s(k,l) = median(tmp);
        end
        
    end
    bold12 = median(data{elec}.allbootsS12,2);
    bold34 = median(data{elec}.allbootsS34,2);
    r2_data_data(l) = ns_cod(bold34,bold12, true);

end

%% V1: now plot simulation LFP and BOLD output versus data for all electrodes

cm = lines(length(find(ismember(v_area,1))));

figure('Position',[0 0 700 400])
for k = 1:8
    subplot(2,4,k),hold on
    signal_use = 4; % bold
    x = squeeze(all_simulation(signal_use,ismember(v_area,[1]),:,k));% signal, all electrodes all conditions
    y = squeeze(all_data(signal_use,ismember(v_area,[1]),:)); % signal, all electrodes all conditions
    for l = 1:size(x,1)% regression line for each electrode
        x_el = x(l,:); y_el = y(l,:);
        p = polyfit(x_el(~isnan(x_el)),y_el(~isnan(x_el)),1);
        x_line=[min(x(:)):0.001:max(x(:))];
        plot(x_line,p(1)*x_line + p(2),'Color',cm(l,:))
        plot(x_el,y_el,'.','MarkerSize',10,'Color',cm(l,:))
        clear x_el y_el % housekeeping
    end
    
    xlim([min(x(:)) max(x(:))]),ylim([min(y(:)) max(y(:))])
    title(['mean R^2 = ' num2str(median(r2_data_fit(k,ismember(v_area,[1]))),2)])
    xlabel('simulated bold'),ylabel('measured bold')   
end
set(gcf,'PaperPositionMode','auto')
print('-depsc','-r300',fullfile(boldlfp_RootPath,'figures', 'FigS7_V1'))
print('-dpng','-r300',fullfile(boldlfp_RootPath,'figures', 'FigS7_V1'))

figure('Position',[0 0 700 400])
ecog_colors = {'k','m','g'};
for k = 1:8
    subplot(2,4,k),hold on
    plot([-1 1],[-1 1],'k')
    for signal_use = 1:3
        x = squeeze(all_simulation(signal_use,ismember(v_area,[1]),:,k));% bb, all electrodes all conditions
        y = squeeze(all_data(signal_use,ismember(v_area,[1]),:)); % bb, all electrodes all conditions
        plot(x,y,'.','MarkerSize',10,'Color',ecog_colors{signal_use})
    end
    axis tight
    xlabel('simulated lfp'),ylabel('measured lfp')   
end

set(gcf,'PaperPositionMode','auto')
print('-depsc','-r300',fullfile(boldlfp_RootPath,'figures', 'FigS9A_V1'))
print('-dpng','-r300',fullfile(boldlfp_RootPath,'figures', 'FigS9A_V1'))


%% V2/V3: now plot simulation LFP and BOLD output versus data for all electrodes
cm = lines(length(find(ismember(v_area,[2 3]))));

figure('Position',[0 0 700 400])
for k = 1:8
    subplot(2,4,k),hold on
    signal_use = 4; % bold
    x = squeeze(all_simulation(signal_use,ismember(v_area,[2 3]),:,k));% signal, all electrodes all conditions
    y = squeeze(all_data(signal_use,ismember(v_area,[2 3]),:)); % signal, all electrodes all conditions
    for l = 1:size(x,1)% regression line for each electrode
        x_el = x(l,:); y_el = y(l,:);
        p = polyfit(x_el(~isnan(x_el)),y_el(~isnan(x_el)),1);
        x_line=[min(x(:)):0.001:max(x(:))];
        plot(x_line,p(1)*x_line + p(2),'Color',cm(l,:))
        plot(x_el,y_el,'.','MarkerSize',10,'Color',cm(l,:))
        clear x_el y_el % housekeeping
    end
    xlim([min(x(:)) max(x(:))]),ylim([min(y(:)) max(y(:))])
    title(['Median R^2 = ' num2str(median(r2_data_fit(k,ismember(v_area,[2 3]))),2)])
    xlabel('simulated bold'),ylabel('measured bold')   
end
set(gcf,'PaperPositionMode','auto')
print('-depsc','-r300',fullfile(boldlfp_RootPath,'figures', 'FigS7_V2V3'))
print('-dpng','-r300',fullfile(boldlfp_RootPath,'figures', 'FigS7_V2V3'))

ecog_colors = {'k','m','g'};
figure('Position',[0 0 700 400])
for k = 1:8
    subplot(2,4,k),hold on
    plot([-1 1],[-1 1],'k')
    for signal_use = 1:3
        x = squeeze(all_simulation(signal_use,ismember(v_area,[2 3]),:,k));% bb, all electrodes all conditions
        y = squeeze(all_data(signal_use,ismember(v_area,[2 3]),:)); % bb, all electrodes all conditions
        plot(x,y,'.','MarkerSize',10,'Color',ecog_colors{signal_use})
    end
    axis tight
    xlabel('simulated lfp'),ylabel('measured lfp')   
end

set(gcf,'PaperPositionMode','auto')
print('-depsc','-r300',fullfile(boldlfp_RootPath,'figures', 'FigS9B_V2V3'))
print('-dpng','-r300',fullfile(boldlfp_RootPath,'figures', 'FigS9B_V2V3'))


%% Make Figure 7C: for one model all r2 simulated versus measured BOLD

figure('Position',[0 0 300 200]),hold on
model_plot = 1;

subplot(1,2,1),hold on
y = r2_data_fit(model_plot,v_area==1);
boxplot(y,'Width',.4);
hold on

% plot R2 from reshuffeling
plot([1-.4 1+.4],median(r2_data_fit_s(model_plot,v_area==1)) * [1 1], ...
    ':','Color',[.5 .5 .5],'LineWidth',2)

% plot R2 from test-retest
plot([1-.4 1+.4],median(r2_data_data(v_area==1)) * [1 1], ...
 '-','Color',[.5 .5 .5],'LineWidth',2)

plot([1-.1:.2./(length(y)-1):1+.1],y,'k.','MarkerSize',20)
plot([1-.1:.2./(length(y)-1):1+.1],y,'y.','MarkerSize',10)

ylim([-1.1 1])
set(gca, 'YTick',-1:.5:1);
ylabel('R^2')

subplot(1,2,2),hold on
y = r2_data_fit(model_plot,v_area==2 | v_area==3);
boxplot(y,'Width',.4);
plot([1-.1:.2./(length(y)-1):1+.1],y,'k.','MarkerSize',20)
plot([1-.1:.2./(length(y)-1):1+.1],y,'y.','MarkerSize',10)
ylim([-1.1 1])

hold on
% plot R2 from reshuffeling
plot([1-.4 1+.4],median(r2_data_fit_s(model_plot,v_area==2 | v_area==3)) * [1 1], ...
    ':','Color',[.5 .5 .5],'LineWidth',2)

% plot R2 from test-retest
plot([1-.4 1+.4],median(r2_data_data(v_area==2 | v_area==3)) * [1 1], ...
 '-','Color',[.5 .5 .5],'LineWidth',2)
ylabel('R^2')

set(gcf,'PaperPositionMode','auto')
print('-depsc','-r300',fullfile(boldlfp_RootPath,'figures', 'Fig7C'))
print('-dpng','-r300',fullfile(boldlfp_RootPath,'figures', 'Fig7C'))


