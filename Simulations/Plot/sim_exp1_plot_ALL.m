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
data_DTI_WM1 = importdata(['/Users/live/Documents/Masterthesis/MATLAB/DKI/Results_corrected/simulation_experiment_1/NSA=' num2str(NSA) '/Pool/voxel1_DTI.txt']);
fa_DTI_WM1 = data_DTI_WM1(1:3000,:);
md_DTI_WM1 = data_DTI_WM1(3001:6000,:);
rd_DTI_WM1 = data_DTI_WM1(6001:9000,:);
ad_DTI_WM1 = data_DTI_WM1(9001:12000,:);
% DTI analytic
data_DTI2_WM = importdata(['/Users/live/Documents/Masterthesis/MATLAB/DKI/Results_corrected/simulation_experiment_1/NSA=' num2str(NSA) '/Pool/subject002_voxel1_DTI2_rel.txt']);
fa_DTI2_WM = data_DTI2_WM(1:1000,:);  fa_DTI2_WM(1001:3000,:) = NaN;
md_DTI2_WM = data_DTI2_WM(1001:2000,:); md_DTI2_WM(1001:3000,:) = NaN;
rd_DTI2_WM = data_DTI2_WM(2001:3000,:); rd_DTI2_WM(1001:3000,:) = NaN;
ad_DTI2_WM = data_DTI2_WM(3001:4000,:); ad_DTI2_WM(1001:3000,:) = NaN;


fa_pool = [fa_DKI_WM fa_DTI_WM1 fa_DTI2_WM];
md_pool = [md_DKI_WM md_DTI_WM1 md_DTI2_WM];
rd_pool = [rd_DKI_WM rd_DTI_WM1 rd_DTI2_WM];
ad_pool = [ad_DKI_WM ad_DTI_WM1 ad_DTI2_WM];
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
% DTI, WLS b 1000
data_DTI_GM1 = importdata(['/Users/live/Documents/Masterthesis/MATLAB/DKI/Results_corrected/simulation_experiment_1/NSA=' num2str(NSA) '/Pool/voxel2_DTI.txt']);
fa_DTI_GM1 = data_DTI_GM1(1:3000,:);
md_DTI_GM1 = data_DTI_GM1(3001:6000,:);
rd_DTI_GM1 = data_DTI_GM1(6001:9000,:);
ad_DTI_GM1 = data_DTI_GM1(9001:12000,:);
% DTI analytic
data_DTI_GM3 = importdata(['/Users/live/Documents/Masterthesis/MATLAB/DKI/Results_corrected/simulation_experiment_1/NSA=' num2str(NSA) '/Pool/subject002_voxel2_DTI2_rel.txt']);
fa_DTI_GM3 = data_DTI_GM3(1:1000,:);  fa_DTI_GM3(1001:3000,:) = NaN;
md_DTI_GM3 = data_DTI_GM3(1001:2000,:); md_DTI_GM3(1001:3000,:) = NaN;
rd_DTI_GM3 = data_DTI_GM3(2001:3000,:); rd_DTI_GM3(1001:3000,:) = NaN;
ad_DTI_GM3 = data_DTI_GM3(3001:4000,:); ad_DTI_GM3(1001:3000,:) = NaN;


fa_pool = [fa_DKI_GM fa_DTI_GM1 fa_DTI_GM3];
md_pool = [md_DKI_GM md_DTI_GM1 md_DTI_GM3];
rd_pool = [rd_DKI_GM rd_DTI_GM1 rd_DTI_GM3];
ad_pool = [ad_DKI_GM ad_DTI_GM1 ad_DTI_GM3];
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
% DTI, WLS b 1000 
data_DTI_CSF1 = importdata(['/Users/live/Documents/Masterthesis/MATLAB/DKI/Results_corrected/simulation_experiment_1/NSA=' num2str(NSA) '/Pool/voxel3_DTI.txt']);
fa_DTI_CSF1 = data_DTI_CSF1(1:3000,:);
md_DTI_CSF1 = data_DTI_CSF1(3001:6000,:);
rd_DTI_CSF1 = data_DTI_CSF1(6001:9000,:);
ad_DTI_CSF1 = data_DTI_CSF1(9001:12000,:);
% DTI analytic
data_DTI_CSF3 = importdata(['/Users/live/Documents/Masterthesis/MATLAB/DKI/Results_corrected/simulation_experiment_1/NSA=' num2str(NSA) '/Pool/subject002_voxel3_DTI2_rel.txt']);
fa_DTI_CSF3 = data_DTI_CSF3(1:1000,:);  fa_DTI_CSF3(1001:3000,:) = NaN;
md_DTI_CSF3 = data_DTI_CSF3(1001:2000,:); md_DTI_CSF3(1001:3000,:) = NaN;
rd_DTI_CSF3 = data_DTI_CSF3(2001:3000,:); rd_DTI_CSF3(1001:3000,:) = NaN;
ad_DTI_CSF3 = data_DTI_CSF3(3001:4000,:); ad_DTI_CSF3(1001:3000,:) = NaN;


fa_pool = [fa_DKI_CSF fa_DTI_CSF1 fa_DTI_CSF3];
md_pool = [md_DKI_CSF md_DTI_CSF1 md_DTI_CSF3];
rd_pool = [rd_DKI_CSF rd_DTI_CSF1 rd_DTI_CSF3];
ad_pool = [ad_DKI_CSF ad_DTI_CSF1 ad_DTI_CSF3];
mk_pool = [mk_DKI_CSF];
rk_pool = [rk_DKI_CSF];
ak_pool = [ak_DKI_CSF];


%% Plot DKI and DTI1 and DTI2 together (groups of 3)

params_pool = [{fa_pool} {md_pool} {rd_pool} {ad_pool}]; % DKI DTI DTI2
params_string = ['FA'; 'MD'; 'RD'; 'AD'];

group = [1 4 7 10 13 16 19 22 25 28 2 5 8 11 14 17 20 23 26 29 3 6 9 12 15 18 21 24 27 30];
positions = [1 1.25 1.5 2.25 2.5 2.75 3.5 3.75 4 4.75 5 5.25 6 6.25 6.5 7.25 7.5 7.75 8.5 8.75 9 9.75 10 10.25 11 11.25 11.5 ...
    12.25 12.5 12.75];
folder_pool = ['/Users/live/Documents/Masterthesis/MATLAB/DKI/Results_corrected/simulation_experiment_1/NSA=' num2str(NSA) '/Pool/Plots_together'];
fontsize = 20;

color1 = [0, 0.4470, 0.7410]; % DKI color (blue)
color2 = [0.8500, 0.3250, 0.0980]; % DTI1 color (orange)
color3 = [0.55, 0.71, 0]; % DTI2 color (green)

color = [{color3},{color2},{color1},{color3},{color2},{color1},{color3},{color2},{color1},{color3},{color2},{color1},{color3},{color2},{color1},...
    {color3},{color2},{color1},{color3},{color2},{color1},{color3},{color2},{color1},{color3},{color2},{color1},{color3},{color2},{color1}];

% x axis
snr = {'10','15','20','25','30','40','50','60','80','100'};

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
    
    % Legend
    set(gca, 'fontsize', fontsize)
    c = get(gca, 'Children');
    hleg1 = legend(c(1:3), 'DKI', 'DTI', 'DTI, analytic');
    %legend('location','southeast');


    xlabel('SNR', 'fontsize', fontsize)
    ylabel('Relative error [%]', 'fontsize', fontsize);
    title([params_string(i,:)], 'fontsize', fontsize)
    
    %pause;
    w = 12;
    h = 6;
    set(gcf, 'PaperSize', [w h]);
    set(gcf, 'PaperPosition', [0 0 w h]);
    saveas(figure(i), fullfile(folder_pool, ['simexp1_' tissue '_DTI2_' params_string(i,:)]), 'pdf');

end


