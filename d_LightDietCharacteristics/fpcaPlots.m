%% Plot fPCA 
% Lara Weed
% 19 May 2025

%% Load Data

ft = readtable('/Users/lara/Library/CloudStorage/OneDrive-Stanford/Research/Projects/Light Circadian/Data/FPCA/fpca_mean_and_sd_components_light_lara.csv');

vn = ft.Properties.VariableNames;

figure('Renderer','painters','Position',[500 500 700 600])
for i = 1:4
    ax(i) = subplot(2,2,i);
    plot(ft.hour,ft.mean,'k','linewidth',2)
    hold on
    plot(ft.hour,ft.(vn{4+2*(i-1)}),'r','linewidth',2)
    plot(ft.hour,ft.(vn{5+2*(i-1)}),'b','linewidth',2)
    grid on
    set(gca,'fontweight','bold','fontsize',12)
    xticks([0:4:24])
end

linkaxes(ax,'xy')
xlim([0 24])
ylim([-100 700])



