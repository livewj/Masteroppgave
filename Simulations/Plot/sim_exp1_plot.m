% Plots results of simulation experiment 1 as boxplots

%% Load data: White matter
NSA = 1;
tissue = 'WM';

% DKI, WLS
data_DKI_WM = importdata(['/Users/live/Documents/Masterthesis/MATLAB/DKI/Results_corrected/simulation_experiment_1/NSA=' num2str(NSA) '/Pool/voxel1_DKI.txt']);
fa_DKI_WM = data_DKI_WM(1:3000,:);
md_DKI_WM = data_DKI_WM(3001:6000,:);
rd_DKI_WM = data_DKI_WM(6001:9000,:);
ad_DKI_WM = data_DKI_WM(9001:12000,:);
mk_DKI_WM = data_DKI_WM(12001:15000,:);
rk_DKI_WM = data_DKI_WM(15001:18000,:);
ak_DKI_WM = data_DKI_WM(18001:21000,:);
% DTI, WLS b 1000
data_DTI_WM = importdata(['/Users/live/Documents/Masterthesis/MATLAB/DKI/Results_corrected/simulation_experiment_1/NSA=' num2str(NSA) '/Pool/voxel1_DTI.txt']);
fa_DTI_WM = data_DTI_WM(1:3000,:);
md_DTI_WM = data_DTI_WM(3001:6000,:);
rd_DTI_WM = data_DTI_WM(6001:9000,:);
ad_DTI_WM = data_DTI_WM(9001:12000,:);


fa_pool = [fa_DKI_WM fa_DTI_WM];
md_pool = [md_DKI_WM md_DTI_WM];
rd_pool = [rd_DKI_WM rd_DTI_WM];
ad_pool = [ad_DKI_WM ad_DTI_WM];
mk_pool = [mk_DKI_WM];
rk_pool = [rk_DKI_WM];
ak_pool = [ak_DKI_WM];


%% Load data: Gray matter
NSA = 1;
tissue = 'GM';

% DKI, WLS
data_DKI_GM = importdata(['/Users/live/Documents/Masterthesis/MATLAB/DKI/Results_corrected/simulation_experiment_1/NSA=' num2str(NSA) '/Pool/voxel2_DKI.txt']);
fa_DKI_GM = data_DKI_GM(1:3000,:);
md_DKI_GM = data_DKI_GM(3001:6000,:);
rd_DKI_GM = data_DKI_GM(6001:9000,:);
ad_DKI_GM = data_DKI_GM(9001:12000,:);
mk_DKI_GM = data_DKI_GM(12001:15000,:);
rk_DKI_GM = data_DKI_GM(15001:18000,:);
ak_DKI_GM = data_DKI_GM(18001:21000,:);
% DTI, WLS
data_DTI_GM = importdata(['/Users/live/Documents/Masterthesis/MATLAB/DKI/Results_corrected/simulation_experiment_1/NSA=' num2str(NSA) '/Pool/voxel2_DTI.txt']);
fa_DTI_GM = data_DTI_GM(1:3000,:);
md_DTI_GM = data_DTI_GM(3001:6000,:);
rd_DTI_GM = data_DTI_GM(6001:9000,:);
ad_DTI_GM = data_DTI_GM(9001:12000,:);


fa_pool = [fa_DKI_GM fa_DTI_GM];
md_pool = [md_DKI_GM md_DTI_GM];
rd_pool = [rd_DKI_GM rd_DTI_GM];
ad_pool = [ad_DKI_GM ad_DTI_GM];
mk_pool = [mk_DKI_GM];
rk_pool = [rk_DKI_GM];
ak_pool = [ak_DKI_GM];


%% Load data: CSF
NSA = 1;
tissue = 'CSF';

% DKI, WLS
data_DKI_CSF = importdata(['/Users/live/Documents/Masterthesis/MATLAB/DKI/Results_corrected/simulation_experiment_1/NSA=' num2str(NSA) '/Pool/voxel3_DKI.txt']);
fa_DKI_CSF = data_DKI_CSF(1:3000,:);
md_DKI_CSF = data_DKI_CSF(3001:6000,:);
rd_DKI_CSF = data_DKI_CSF(6001:9000,:);
ad_DKI_CSF = data_DKI_CSF(9001:12000,:);
mk_DKI_CSF = data_DKI_CSF(12001:15000,:);
rk_DKI_CSF = data_DKI_CSF(15001:18000,:);
ak_DKI_CSF = data_DKI_CSF(18001:21000,:);
% DTI, WLS
data_DTI_CSF = importdata(['/Users/live/Documents/Masterthesis/MATLAB/DKI/Results_corrected/simulation_experiment_1/NSA=' num2str(NSA) '/Pool/voxel3_DTI.txt']);
fa_DTI_CSF = data_DTI_CSF(1:3000,:);
md_DTI_CSF = data_DTI_CSF(3001:6000,:);
rd_DTI_CSF = data_DTI_CSF(6001:9000,:);
ad_DTI_CSF = data_DTI_CSF(9001:12000,:);


fa_pool = [fa_DKI_CSF fa_DTI_CSF];
md_pool = [md_DKI_CSF md_DTI_CSF];
rd_pool = [rd_DKI_CSF rd_DTI_CSF];
ad_pool = [ad_DKI_CSF ad_DTI_CSF];
mk_pool = [mk_DKI_CSF];
rk_pool = [rk_DKI_CSF];
ak_pool = [ak_DKI_CSF];


%% Plot DKI and DTI together

params_pool = [{fa_pool} {md_pool} {rd_pool} {ad_pool}];
params_string = ['FA'; 'MD'; 'RD'; 'AD'];
group = [1 3 5 7 9 11 13 15 17 19 2 4 6 8 10 12 14 16 18 20];
positions = [1 1.25 2 2.25 3 3.25 4 4.25 5 5.25 6 6.25 7 7.25 8 8.25 9 9.25 10 10.25];
folder_pool = ['/Users/live/Documents/Masterthesis/MATLAB/DKI/Results_corrected/simulation_experiment_1/NSA=' num2str(NSA) '/Pool/Plots_together'];
fontsize = 20;

color1 = [0.8500, 0.3250, 0.0980]; % DTI color (orange)
color2 = [0, 0.4470, 0.7410]; % DKI color (blue)
color = [{color1},{color2},{color1},{color2},{color1},{color2},{color1},{color2},{color1},{color2},{color1},{color2},{color1},{color2},{color1},{color2},{color1},{color2},{color1},{color2}];
    
% x axis
snr = {'10','15','20','25','30','40','50','60','80','100'};

for i = 1:length(params_pool)
    x = params_pool{i};
    figure(i), boxplot(x, group, 'positions', positions);
    
    hold on
    plot(linspace(0,105), 0*ones(1,100), '-');
    hold off
    
    % Group name
    set(gca,'xtick',[mean(positions(1:2)) mean(positions(3:4)) mean(positions(5:6)) mean(positions(7:8)) ...
        mean(positions(9:10)) mean(positions(11:12)) mean(positions(13:14)) mean(positions(15:16)) ...
        mean(positions(17:18)) mean(positions(19:20))])
    set(gca,'xticklabel',snr)
    
    % Change box colors
    h = findobj(gca,'Tag','Box');
    for j=1:length(h)
       patch(get(h(j),'XData'),get(h(j),'YData'),color{j},'FaceAlpha',.5);
    end
    
    % Change color of outliers
    outliers = findobj(gcf,'tag', 'Outliers');
    for j = 1:length(outliers)
        outliers(j).MarkerEdgeColor = color{j};
    end
    
    % Median color
    % lines = findobj(gcf, 'type', 'line', 'Tag', 'Median');
    % set(lines, 'Color', 'g');
    
    % Legend
    set(gca, 'fontsize', fontsize)
    c = get(gca, 'Children');
    hleg1 = legend(c(1:2), 'DKI', 'DTI');
    %legend('location','southeast');


    xlabel('SNR', 'fontsize', fontsize)
    ylabel('Relative error [%]', 'fontsize', fontsize);
    title([params_string(i,:)], 'fontsize', fontsize)
    
    %pause;
    w = 12;
    h = 6;
    set(gcf, 'PaperSize', [w h]);
    set(gcf, 'PaperPosition', [0 0 w h]);
    saveas(figure(i), fullfile(folder_pool, ['simexp1_' tissue '_' params_string(i,:)]), 'pdf');

end

%% Plot only the DKI parameters

params_pool = [{mk_pool} {rk_pool} {ak_pool}];
params_string = ['MK'; 'RK'; 'AK'];
folder_pool = ['/Users/live/Documents/Masterthesis/MATLAB/DKI/Results_corrected/simulation_experiment_1/NSA=' num2str(NSA) '/Pool/Plots_together'];
fontsize = 20;


for i = 1:length(params_pool)
    figure(i), clf, boxplot(params_pool{i}, SNR)
    hold on
    plot(linspace(0,105), 0*ones(1,100), '-');
    hold off
    
    % Change box colors
    h = findobj(gca,'Tag','Box');
    for j=1:length(h)
       patch(get(h(j),'XData'),get(h(j),'YData'),color2,'FaceAlpha',.5);
    end
    
    % Change color of outliers
    outliers = findobj(gcf,'tag', 'Outliers');
    for j = 1:length(outliers)
        outliers(j).MarkerEdgeColor = color2;
    end
    
    % Legend
    set(gca, 'fontsize', fontsize)
    c = get(gca, 'Children');
    hleg1 = legend(c(1:1), 'DKI'); % XXX endre
    
    xlabel('SNR', 'fontsize', fontsize)
    ylabel('Relative error [%]', 'fontsize', fontsize);
    title([params_string(i,:)], 'fontsize', fontsize)
    
    %pause;
    w = 7;
    h = 5;
    set(gcf, 'PaperSize', [w h]);
    set(gcf, 'PaperPosition', [0 0 w h]);
    saveas(figure(i), fullfile(folder_pool, ['simexp1_' tissue '_' params_string(i,:)]), 'pdf');
    
end

%% Plot MK, RK, AK for WM, GM and CSF in SAME PLOT (groups of 3)

rk_DKI_WM(2001:end,:) = NaN;

mk_pool = [mk_DKI_WM mk_DKI_GM mk_DKI_CSF];
rk_pool = [rk_DKI_WM rk_DKI_GM rk_DKI_CSF];
ak_pool = [ak_DKI_WM ak_DKI_GM ak_DKI_CSF];

params_pool = [{mk_pool} {rk_pool} {ak_pool}]; % DKI
folder_pool = ['/Users/live/Documents/Masterthesis/MATLAB/DKI/Results_corrected/simulation_experiment_1/NSA=1/Pool/Plots_together'];
params_string = ['MK'; 'RK'; 'AK'];
group = [1 4 7 10 13 16 19 22 25 28 2 5 8 11 14 17 20 23 26 29 3 6 9 12 15 18 21 24 27 30];
positions = [1 1.25 1.5 2.25 2.5 2.75 3.5 3.75 4 4.75 5 5.25 6 6.25 6.5 7.25 7.5 7.75 8.5 8.75 9 9.75 10 10.25 11 11.25 11.5 ...
    12.25 12.5 12.75];

colorWM =  [0.98, 0.91, 0.71]; % white
colorGM = [0.41, 0.41, 0.41]; % gray
colorCSF = [0.34, 0.63, 0.83]; % blue

color = [{colorCSF},{colorGM},{colorWM},{colorCSF},{colorGM},{colorWM},{colorCSF},{colorGM},{colorWM},{colorCSF},{colorGM},{colorWM},{colorCSF},{colorGM},{colorWM},...
    {colorCSF},{colorGM},{colorWM},{colorCSF},{colorGM},{colorWM},{colorCSF},{colorGM},{colorWM},{colorCSF},{colorGM},{colorWM},{colorCSF},{colorGM},{colorWM}];

for i = 1:length(params_pool)
    x = params_pool{i};
    figure(i), boxplot(x, group, 'positions', positions);
    
    hold on
    plot(linspace(0,105), 0*ones(1,100), '-');
    hold off

    % Group name
    set(gca,'xtick',[mean(positions(1:3)) mean(positions(4:6)) mean(positions(7:9)) mean(positions(10:12)) ...
        mean(positions(13:15)) mean(positions(16:18)) mean(positions(19:21)) mean(positions(22:24)) ...
        mean(positions(25:27)) mean(positions(28:30))])
    set(gca,'xticklabel',{'10','15','20','25','30','40','50','60','80','100'})


    % Change box colors
    h = findobj(gca,'Tag','Box');
    for j=1:length(h)
       patch(get(h(j),'XData'),get(h(j),'YData'),color{j},'FaceAlpha',.5);
    end
    
    % Change color of outliers
    outliers = findobj(gcf,'tag', 'Outliers');
    for j = 1:length(outliers)
        outliers(j).MarkerEdgeColor = color{j};
    end
    
    % Legend
    set(gca, 'fontsize', fontsize)
    c = get(gca, 'Children');
    hleg1 = legend(c(1:3),'White matter','Gray matter','CSF');
    %legend('location','southeast');

    xlabel('SNR', 'fontsize', fontsize)
    ylabel('Relative error [%]', 'fontsize', fontsize);
    title([params_string(i,:)], 'fontsize', fontsize)
    
    %pause;
    w = 12;
    h = 6;
    set(gcf, 'PaperSize', [w h]);
    set(gcf, 'PaperPosition', [0 0 w h]);
    saveas(figure(i), fullfile(folder_pool, ['simexp1_DKI_' params_string(i,:)]), 'pdf');

end

%% Plot MD for WM, GM and CSF for DTI and DKI in the SAME PLOT (groups of 6)

md_pool = [md_DKI_WM md_DTI_WM md_DKI_GM md_DTI_GM md_DKI_CSF md_DTI_CSF];

params_pool = [{md_pool}]; % DKI DTI
folder_pool = ['/Users/live/Documents/Masterthesis/MATLAB/DKI/Results_corrected/simulation_experiment_1/NSA=1/Pool/Plots_together'];
params_string = ['MD'];

group = [1 7 13 19 25 31 37 43 49 55 2 8 14 20 26 32 38 44 50 56 3 9 15 21 27 33 39 45 51 57 4 10 16 22 28 34 40 46 52 58 ...
    5 11 17 23 29 35 41 47 53 59 6 12 18 24 30 36 42 48 54 60];
positions = [1 1.25 1.5 1.75 2 2.25 3 3.25 3.5 3.75 4 4.25 5 5.25 5.5 5.75 6 6.25 7 7.25 7.5 7.75 8 8.25 ...
    9 9.25 9.5 9.75 10 10.25 11 11.25 11.5 11.75 12 12.25 13 13.25 13.5 13.75 14 14.25 15 15.25 15.5 15.75 16 16.25...
    17 17.25 17.5 17.75 18 18.25 19 19.25 19.50 19.75 20 20.25];

% DKI
colorWM1 =  [0.98, 0.91, 0.71]; % banana mania
colorGM1 = [0.41, 0.41, 0.41];  % dim gray
colorCSF1 = [0.34, 0.63, 0.83]; % blue

% DTI
colorWM2 =  [1, 0.60, 0.40]; % desert
colorGM2 = [0, 0, 0];        % black
colorCSF2 = [0, 0, 0.9];     % french gray

color = [{colorCSF2},{colorCSF1},{colorGM2},{colorGM1},{colorWM2},{colorWM1},{colorCSF2},{colorCSF1},{colorGM2},{colorGM1},{colorWM2},{colorWM1},{colorCSF2},{colorCSF1},{colorGM2},{colorGM1},{colorWM2},{colorWM1},...
    {colorCSF2},{colorCSF1},{colorGM2},{colorGM1},{colorWM2},{colorWM1},{colorCSF2},{colorCSF1},{colorGM2},{colorGM1},{colorWM2},{colorWM1},{colorCSF2},{colorCSF1},{colorGM2},{colorGM1},{colorWM2},{colorWM1},...
    {colorCSF2},{colorCSF1},{colorGM2},{colorGM1},{colorWM2},{colorWM1},{colorCSF2},{colorCSF1},{colorGM2},{colorGM1},{colorWM2},{colorWM1},{colorCSF2},{colorCSF1},{colorGM2},{colorGM1},{colorWM2},{colorWM1},{colorCSF2},{colorCSF1},{colorGM2},{colorGM1},{colorWM2},{colorWM1}];

for i = 1:length(params_pool)
    x = params_pool{i};
    figure(i), boxplot(x, group, 'positions', positions);
    
    hold on
    plot(linspace(0,105), 0*ones(1,100), '-');
    hold off

    % Group name
    set(gca,'xtick',[mean(positions(1:6)) mean(positions(7:12)) mean(positions(13:18)) mean(positions(19:24)) ...
        mean(positions(25:30)) mean(positions(31:36)) mean(positions(37:42)) mean(positions(43:48)) ...
        mean(positions(49:54)) mean(positions(55:60))])
    set(gca,'xticklabel',{'10','15','20','25','30','40','50','60','80','100'})


    % Change box colors
    h = findobj(gca,'Tag','Box');
    for j=1:length(h)
       patch(get(h(j),'XData'),get(h(j),'YData'),color{j},'FaceAlpha',.5);
    end
    
    % Change color of outliers
    outliers = findobj(gcf,'tag', 'Outliers');
    for j = 1:length(outliers)
        outliers(j).MarkerEdgeColor = color{j};
    end
    
    % Legend
    set(gca, 'fontsize', fontsize)
    c = get(gca, 'Children');
    hleg1 = legend(c(1:6),'White matter, DKI','White matter, DTI', 'Gray matter, DKI', 'Gray matter, DTI',...
        'CSF, DKI', 'CSF, DTI');
    legend('location','southeast');

    xlabel('SNR', 'fontsize', fontsize)
    ylabel('Relative error [%]', 'fontsize', fontsize);
    title([params_string(i,:)], 'fontsize', fontsize)
    
    %pause;
    w = 12;
    h = 6;
    set(gcf, 'PaperSize', [w h]);
    set(gcf, 'PaperPosition', [0 0 w h]);
    saveas(figure(i), fullfile(folder_pool, ['simexp1_DKIDTI_all_' params_string(i,:)]), 'pdf');

end


