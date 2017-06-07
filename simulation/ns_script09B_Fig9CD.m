% Reproduce pannels from Figure 9 from the following paper:
%  
%   Hermes, Nguyen and Winawer (2017). Neuronal synchrony and the relation
%   between the BOLD signal and the local field potential. PLOS Biology
%   http://dx.doi.org/...
%
% Run:
% ns_script09B_Fig9CD
%
% Purpose: load simulated neural data and show the relation between LFP and
% BOLD
%
% DH and JW 2016

%% Load all fitted electrodes

clear all
sim_nr = 2;
els = 1:22;
nr_elec = length(els);

%%% OUTPUTS:
v_area = NaN(length(els),1); %visual area per electrode

% SIMULATION: 4:bb,g,a,bold, nr els, up to 10 conditions, 8 simulations
all_simulation = NaN(4,length(els),10,8); % BOLD simulation

% SIMULATION output regression models
cod_crossval_out = NaN(length(els),7); % r2 for regression models
cod_crossval_outShuff = NaN(length(els),7); % r2 for regression models
all_regressbeta   = NaN(length(els),7,4); % betas for regression models

load(fullfile(boldlfp_RootPath, 'data', 'boldecog_structure_final'));

for l = 1:nr_elec
    
    elec = els(l);
    
    disp(['el ' int2str(l) ' of ' int2str(nr_elec)])
    v_area(l) = data{l}.v_area;
    
    % load the simulation outputs 
    load(fullfile(boldlfp_RootPath, 'data', sprintf('NS_simnr%d_elec%d_simulation_outputs',sim_nr,elec)),'simulation_outputs')
       
    % get simulated ECoG (bb, g, a) and BOLD responses into
    % 'all_simulation' matrix
    all_simulation(:,l,1:size(simulation_outputs,1),:) = permute(simulation_outputs, [3 1 2]);
    
    % load output from the first model (BB - level, G - coh, A - level)
    prm_set = 1;
    load(fullfile(boldlfp_RootPath, 'data', sprintf('NS_simnr%d_elec%d_NS_prmset%d',sim_nr,elec,prm_set)),'NS')
    
    % make sure we are using COD rather than r2 for accuracy measure
    NS = ns_summary_statistics(NS); %disp(NS.stats) 
    
    for k = 1:length(NS.stats)
        % cross validated R2:
        cod_crossval_out(l,k) = median(NS.stats(k).stats(:,3));
        cod_crossval_outShuff(l,k) = median(NS.stats(k).stats_shuffled(:,3));
        % beta values:
        temp_beta = median(NS.stats(k).beta(:,2:end),1);
        all_regressbeta(l,k,1:length(temp_beta)) = temp_beta;
    end
end


%% R2 plots averaged for V1 and V2 simulations

bar_colors={[1 0 0],[1 1 0],[1 .5 0],[0 .2 1],[.5 0 1],[0 .5 0],[.4 .2 .1]};
box_colors = zeros(length(bar_colors),3);
for ii = 1:length(bar_colors), box_colors(ii,:) = bar_colors{ii}; end

figure('Position',[0 0 580 200])
% ----------------


% CROSS-VALIDATED R^2 when taking all boots
for whichAreas = 1:2
    
    subplot(1,2,whichAreas),hold on % plot V1
    if whichAreas == 1
        whichElectrodes = v_area==1;
    else
        whichElectrodes = v_area==2 | v_area==3;
    end
    
    boxplot(cod_crossval_out(whichElectrodes==1,1:7),'Colors', box_colors);%, 'PlotStyle', 'compact');
    for k=1:size(cod_crossval_out,2)
        
        % plot R2 from reshuffling
        plot([k-.4 k+.4],[median(median(cod_crossval_outShuff(whichElectrodes,k,:),3),1) ...
            median(median(cod_crossval_outShuff(whichElectrodes,k,:),3),1)],':','Color',[.5 .5 .5],'LineWidth',2)
        
    end
    
    clear mean_resp st_err
    xlim([0 8]),ylim([-.8 1])
    set(gca,'XTick',[1:7],'XTickLabel',{'bb','g','bb_g','a','bb_a','g_a','bb_g_a'})
    set(gca,'YTick',[-1:.2:1])
    
    if whichAreas == 1, title('V1 R^2 cross-val')
    else, title('V2/V3 R^2 cross-val'); end
    
end

set(gcf,'PaperPositionMode','auto')
fname = fullfile(boldlfp_RootPath, 'figures', 'Fig9CD');
print('-dpng','-r300',fname)
print('-depsc','-r300',fname) 

disp(['R^2: ' num2str(median(cod_crossval_out(v_area==1,:),1))]);
disp(['R^2: ' num2str(median(cod_crossval_out(v_area==2 | v_area==3,:),1))]);


%% BETA plots averaged for V1 and V2 simulations

labels_beta={{'bb','',''},{'','g',''},{'bb','g',''},{'','','a'},...
    {'bb','','a'},{'','g','a'},{'bb','g','a'}};
labels_index={1,2,[1 2],3,[1 3],[2 3],[1 2 3]};
bb_g_a_color={[.9 .9 .9],[.6 .6 .6],[.3 .3 .3]};

% plot V1
figure('Position',[0 0 450 100])
for k=1:size(all_regressbeta,2)
    xl_ind=labels_index{k};
    subplot(1,size(all_regressbeta,2)*2,k),hold on 
    
    for m=1:3 % nr of betas
        % take the median across the bootstraps for each electrode
        temp_beta=all_regressbeta(v_area==1,k,m);
        if ~isnan(temp_beta)
            % plot mean across electrodes
            bar(xl_ind(m),mean(temp_beta),.7,'FaceColor',bb_g_a_color{xl_ind(m)})
            % plot 2 x standard error as error bar
            st_err = std(temp_beta)./sqrt(length(temp_beta));
            plot([xl_ind(m) xl_ind(m)],[mean(temp_beta)-st_err mean(temp_beta)+st_err],'k')
            % test for significant difference from zero across electrodes using a t-test
            [~,p]=ttest(temp_beta);
            if p<=0.05
                plot(xl_ind(m),-.2,'r*')
            end
        end
    end
    xlim([.5 3.5]),ylim([-.2 .4])
    set(gca,'XTick',1:3,'XTickLabel',labels_beta{k},'YTick',-0.4:.2:.8,'YTickLabel',[])
end

% plot V2/V3
for k=1:size(all_regressbeta,2)
    xl_ind=labels_index{k};
    subplot(1,size(all_regressbeta,2)*2,size(all_regressbeta,2)+k),hold on 

    for m=1:3 % nr of betas
        % take the median across the bootstraps for each electrode
        temp_beta=all_regressbeta(v_area==2 | v_area==3,k,m);
        if ~isnan(temp_beta)
            % plot mean across electrodes
            bar(xl_ind(m),mean(temp_beta),.7,'FaceColor',bb_g_a_color{xl_ind(m)})
            % plot 2 x standard error as error bar
            st_err = std(temp_beta)./sqrt(length(temp_beta));
            plot([xl_ind(m) xl_ind(m)],[mean(temp_beta)-st_err mean(temp_beta)+st_err],'k')
            % test for significant difference from zero across electrodes using a t-test
            [~,p]=ttest(temp_beta);
            if p<=0.05
                plot(xl_ind(m),-.2,'r*')
            end
        end
    end   
    xlim([.5 3.5]),ylim([-.2 .4])
    set(gca,'XTick',1:3,'XTickLabel',labels_beta{k},'YTick',-0.4:.2:.8,'YTickLabel',[])
end

set(gcf,'PaperPositionMode','auto')
fname = fullfile(boldlfp_RootPath, 'figures', 'Fig9CD_betas');
print('-dpng','-r300',fname)
print('-depsc','-r300',fname) 


