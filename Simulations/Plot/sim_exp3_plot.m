% Plots results of simulation experiment 3 as boxplots

%% Load data: White matter
tissue = 'WM';

% DKI
data_DKI_WM = importdata(['/Users/live/Documents/Masterthesis/MATLAB/DKI/Results_corrected/simulation_experiment_3/Pool/voxel1_DKI.txt']);
fa_DKI_WM = data_DKI_WM(1:3000,:);
md_DKI_WM = data_DKI_WM(3001:6000,:);
rd_DKI_WM = data_DKI_WM(6001:9000,:);
ad_DKI_WM = data_DKI_WM(9001:12000,:);
mk_DKI_WM = data_DKI_WM(12001:15000,:);
rk_DKI_WM = data_DKI_WM(15001:18000,:);
ak_DKI_WM = data_DKI_WM(18001:end,:);

% DTI
data_DTI_WM = importdata(['/Users/live/Documents/Masterthesis/MATLAB/DKI/Results_corrected/simulation_experiment_3/Pool/voxel1_DTI.txt']);
fa_DTI_WM = data_DTI_WM(1:2000,:); fa_DTI_WM(2001:3000,:) = NaN;
md_DTI_WM = data_DTI_WM(2001:4000,:); md_DTI_WM(2001:3000,:) = NaN;
rd_DTI_WM = data_DTI_WM(4001:6000,:); rd_DTI_WM(2001:3000,:) = NaN;
ad_DTI_WM = data_DTI_WM(6001:8000,:); ad_DTI_WM(2001:3000,:) = NaN;

% DTI 3000
data_DTI3_WM = importdata(['/Users/live/Documents/Masterthesis/MATLAB/DKI/Results_corrected/simulation_experiment_3/Pool/voxel1_DTI3.txt']);
fa_DTI3_WM = data_DTI3_WM(1:2000,7:8); fa_DTI3_WM(2001:3000,:) = NaN;
md_DTI3_WM = data_DTI3_WM(2001:4000,7:8); md_DTI3_WM(2001:3000,:) = NaN;
rd_DTI3_WM = data_DTI3_WM(4001:6000,7:8); rd_DTI3_WM(2001:3000,:) = NaN;
ad_DTI3_WM = data_DTI3_WM(6001:8000,7:8); ad_DTI3_WM(2001:3000,:) = NaN;


fa_pool = [fa_DKI_WM fa_DTI_WM];
md_pool = [md_DKI_WM md_DTI_WM];
rd_pool = [rd_DKI_WM rd_DTI_WM];
ad_pool = [ad_DKI_WM ad_DTI_WM];
mk_pool = [mk_DKI_WM];
rk_pool = [rk_DKI_WM];
ak_pool = [ak_DKI_WM];

%% Load data: Gray matter
tissue = 'GM';

% DKI
data_DKI_GM = importdata(['/Users/live/Documents/Masterthesis/MATLAB/DKI/Results_corrected/simulation_experiment_3/Pool/voxel2_DKI.txt']);
fa_DKI_GM = data_DKI_GM(1:3000,:);
md_DKI_GM = data_DKI_GM(3001:6000,:);
rd_DKI_GM = data_DKI_GM(6001:9000,:);
ad_DKI_GM = data_DKI_GM(9001:12000,:);
mk_DKI_GM = data_DKI_GM(12001:15000,:);
rk_DKI_GM = data_DKI_GM(15001:18000,:);
ak_DKI_GM = data_DKI_GM(18001:end,:); ak_DKI_GM(3000,:) = NaN;

% DTI
data_DTI_GM = importdata(['/Users/live/Documents/Masterthesis/MATLAB/DKI/Results_corrected/simulation_experiment_3/Pool/voxel2_DTI.txt']);
fa_DTI_GM = data_DTI_GM(1:2000,:); fa_DTI_GM(2001:3000,:) = NaN;
md_DTI_GM = data_DTI_GM(2001:4000,:); md_DTI_GM(2001:3000,:) = NaN;
rd_DTI_GM = data_DTI_GM(4001:6000,:); rd_DTI_GM(2001:3000,:) = NaN;
ad_DTI_GM = data_DTI_GM(6001:8000,:); ad_DTI_GM(2001:3000,:) = NaN;

% DTI 3000
data_DTI3_GM = importdata(['/Users/live/Documents/Masterthesis/MATLAB/DKI/Results_corrected/simulation_experiment_3/Pool/voxel2_DTI3.txt']);
fa_DTI3_GM = data_DTI3_GM(1:2000,7:8); fa_DTI3_GM(2001:3000,:) = NaN;
md_DTI3_GM = data_DTI3_GM(2001:4000,7:8); md_DTI3_GM(2001:3000,:) = NaN;
rd_DTI3_GM = data_DTI3_GM(4001:6000,7:8); rd_DTI3_GM(2001:3000,:) = NaN;
ad_DTI3_GM = data_DTI3_GM(6001:8000,7:8); ad_DTI3_GM(2001:3000,:) = NaN;


fa_pool = [fa_DKI_GM fa_DTI_GM];
md_pool = [md_DKI_GM md_DTI_GM];
rd_pool = [rd_DKI_GM rd_DTI_GM];
ad_pool = [ad_DKI_GM ad_DTI_GM];
mk_pool = [mk_DKI_GM];
rk_pool = [rk_DKI_GM];
ak_pool = [ak_DKI_GM];


%% Load data: CSF
tissue = 'CSF';

% DKI
data_DKI_CSF = importdata(['/Users/live/Documents/Masterthesis/MATLAB/DKI/Results_corrected/simulation_experiment_3/Pool/voxel3_DKI.txt']);
fa_DKI_CSF = data_DKI_CSF(1:3000,:);
md_DKI_CSF = data_DKI_CSF(3001:6000,:);
rd_DKI_CSF = data_DKI_CSF(6001:9000,:);
ad_DKI_CSF = data_DKI_CSF(9001:12000,:);
mk_DKI_CSF = data_DKI_CSF(12001:15000,:);
rk_DKI_CSF = data_DKI_CSF(15001:18000,:);
ak_DKI_CSF = data_DKI_CSF(18001:end,:);
% DTI
data_DTI_CSF = importdata(['/Users/live/Documents/Masterthesis/MATLAB/DKI/Results_corrected/simulation_experiment_3/Pool/voxel3_DTI.txt']);
fa_DTI_CSF = data_DTI_CSF(1:2000,:); fa_DTI_CSF(2001:3000,:) = NaN;
md_DTI_CSF = data_DTI_CSF(2001:4000,:); md_DTI_CSF(2001:3000,:) = NaN;
rd_DTI_CSF = data_DTI_CSF(4001:6000,:); rd_DTI_CSF(2001:3000,:) = NaN;
ad_DTI_CSF = data_DTI_CSF(6001:8000,:); ad_DTI_CSF(2001:3000,:) = NaN;

% DTI 3000
data_DTI3_CSF = importdata(['/Users/live/Documents/Masterthesis/MATLAB/DKI/Results_corrected/simulation_experiment_3/Pool/voxel3_DTI3.txt']);
fa_DTI3_CSF = data_DTI3_CSF(1:2000,7:8); fa_DTI3_CSF(2001:3000,:) = NaN;
md_DTI3_CSF = data_DTI3_CSF(2001:4000,7:8); md_DTI3_CSF(2001:3000,:) = NaN;
rd_DTI3_CSF = data_DTI3_CSF(4001:6000,7:8); rd_DTI3_CSF(2001:3000,:) = NaN;
ad_DTI3_CSF = data_DTI3_CSF(6001:8000,7:8); ad_DTI3_CSF(2001:3000,:) = NaN;


fa_pool = [fa_DKI_CSF fa_DTI_CSF];
md_pool = [md_DKI_CSF md_DTI_CSF];
rd_pool = [rd_DKI_CSF rd_DTI_CSF];
ad_pool = [ad_DKI_CSF ad_DTI_CSF];
mk_pool = [mk_DKI_CSF];
rk_pool = [rk_DKI_CSF];
ak_pool = [ak_DKI_CSF];

%% One parameter, all sets, DKI

voxel = 'voxel3';

folder_pool = ['/Users/live/Documents/Masterthesis/MATLAB/DKI/Results_corrected/simulation_experiment_3/Pool/'];
fontsize = 20;

model = 'DKI';


params_pool = [{fa_pool(:,1:8)} {md_pool(:,1:8)} {rd_pool(:,1:8)} {ad_pool(:,1:8)} {mk_pool(:,1:8)} {rk_pool(:,1:8)} {ak_pool(:,1:8)}];    
params_string = ['FA'; 'MD'; 'RD'; 'AD'; 'MK'; 'RK'; 'AK'];
color = [0, 0.4470, 0.7410]; % DKI color (blue)

for p = 1:length(params_string)
    one_param_all_sets = [params_pool{p}(:,1) params_pool{p}(:,2) params_pool{p}(:,3) params_pool{p}(:,4) params_pool{p}(:,5) params_pool{p}(:,6) params_pool{p}(:,7) params_pool{p}(:,8)];
    figure(p), clf, boxplot(one_param_all_sets)
    hold on
    plot(linspace(0,105), 0*ones(1,100), '-');
    hold off
    
    % Change box colors
    h = findobj(gca,'Tag','Box');
    for j=1:length(h)
       patch(get(h(j),'XData'),get(h(j),'YData'),color,'FaceAlpha',.5);
    end
    
    % Change color of outliers
    outliers = findobj(gcf,'tag', 'Outliers');
    for j = 1:length(outliers)
        outliers(j).MarkerEdgeColor = color;
    end
    
    % Legend
    set(gca, 'fontsize', fontsize)
    c = get(gca, 'Children');
    hleg1 = legend(c(1:1), model);
    set(gca, 'fontsize', fontsize)
    
    xlabel('Gradient set', 'fontsize', fontsize)
    ylabel('Relative error [%]', 'fontsize', fontsize);
    title([params_string(p,:)], 'fontsize', fontsize)
    
    %pause;
    w = 12;
    h = 6;
    set(gcf, 'PaperSize', [w h]);
    set(gcf, 'PaperPosition', [0 0 w h]);
    saveas(figure(p), fullfile(folder_pool, [voxel '_' model '_' params_string(p,:)]), 'pdf');
end

%% One parameter, all sets, DTI

voxel = 'voxel3';

folder_pool = ['/Users/live/Documents/Masterthesis/MATLAB/DKI/Results_corrected/simulation_experiment_3/Pool/' ];

% fa_pool = [fa_DTI_WM(:,1:6) fa_DTI3_WM];
% md_pool = [md_DTI_WM(:,1:6) md_DTI3_WM];
% rd_pool = [rd_DTI_WM(:,1:6) rd_DTI3_WM];
% ad_pool = [ad_DTI_WM(:,1:6) ad_DTI3_WM];
%  
% fa_pool = [fa_DTI_GM(:,1:6) fa_DTI3_GM];
% md_pool = [md_DTI_GM(:,1:6) md_DTI3_GM];
% rd_pool = [rd_DTI_GM(:,1:6) rd_DTI3_GM];
% ad_pool = [ad_DTI_GM(:,1:6) ad_DTI3_GM];

fa_pool = [fa_DTI_CSF(:,1:6) fa_DTI3_CSF];
md_pool = [md_DTI_CSF(:,1:6) md_DTI3_CSF];
rd_pool = [rd_DTI_CSF(:,1:6) rd_DTI3_CSF];
ad_pool = [ad_DTI_CSF(:,1:6) ad_DTI3_CSF];

model = 'DTI_DTI3';
params_pool = [{fa_pool} {md_pool} {rd_pool} {ad_pool}];    
params_string = ['FA'; 'MD'; 'RD'; 'AD'];
color = [0.8500, 0.3250, 0.0980]; % DTI color (orange) for NSA = 1

for p = 1:length(params_string)
    %one_param_all_sets = [params_pool{p}(:,1) params_pool{p}(:,2) params_pool{p}(:,3) params_pool{p}(:,4) params_pool{p}(:,5) params_pool{p}(:,6) params_pool{p}(:,7) params_pool{p}(:,8)];
    one_param_all_sets = [params_pool{p}(:,1) params_pool{p}(:,2) params_pool{p}(:,3) params_pool{p}(:,4) params_pool{p}(:,5) params_pool{p}(:,6)];
    figure(p), clf, boxplot(one_param_all_sets)
    hold on
    plot(linspace(0,105), 0*ones(1,100), '-');
    hold off
    
    % Change box colors
    h = findobj(gca,'Tag','Box');
    for j=1:length(h)
       patch(get(h(j),'XData'),get(h(j),'YData'),color,'FaceAlpha',.5);
    end
    
    % Change color of outliers
    outliers = findobj(gcf,'tag', 'Outliers');
    for j = 1:length(outliers)
        outliers(j).MarkerEdgeColor = color;
    end
    
    % Legend
    set(gca, 'fontsize', fontsize)
    c = get(gca, 'Children');
    hleg1 = legend(c(1:1), 'DTI');
    set(gca, 'fontsize', fontsize)
    
    xlabel('Gradient set', 'fontsize', fontsize)
    ylabel('Relative error [%]', 'fontsize', fontsize);
    title([params_string(p,:)], 'fontsize', fontsize)
    
    %pause;
    w = 12;
    h = 6;
    set(gcf, 'PaperSize', [w h]);
    set(gcf, 'PaperPosition', [0 0 w h]);
    saveas(figure(p), fullfile(folder_pool, [voxel '_' model '_' params_string(p,:) '_2']), 'pdf');
end

%% MD DTI, all sets, all tissues (groups of 3)

fontsize = 20;
folder_pool = ['/Users/live/Documents/Masterthesis/MATLAB/DKI/Results_corrected/simulation_experiment_3/Pool/'];

%md_pool = [md_DTI_WM md_DTI_GM md_DTI_CSF];
%fa_pool = [fa_DTI_WM fa_DTI_GM fa_DTI_CSF];

%md_pool = [md_DTI_WM(:,1:6) md_DTI3_WM md_DTI_GM(:,1:6) md_DTI3_GM md_DTI_CSF(:,1:6) md_DTI3_CSF];
%pool = [rd_DTI_WM(:,1:6) rd_DTI3_WM rd_DTI_GM(:,1:6) rd_DTI3_GM rd_DTI_CSF(:,1:6) rd_DTI3_CSF];
pool = [ad_DTI_WM(:,1:6) ad_DTI3_WM ad_DTI_GM(:,1:6) ad_DTI3_GM ad_DTI_CSF(:,1:6) ad_DTI3_CSF];
%fa_pool = [fa_DTI_WM(:,1:6) fa_DTI_GM(:,1:6) fa_DTI_CSF(:,1:6)];
   

model = 'DTI';
params_pool = [{pool}];    
params_string = ['AD'];

group = [1 4 7 10 13 16 19 22 2 5 8 11 14 17 20 23 3 6 9 12 15 18 21 24];
positions = [1 1.25 1.5 2.25 2.5 2.75 3.5 3.75 4 4.75 5 5.25 6 6.25 6.5 7.25 7.5 7.75 8.5 8.75 9 9.75 10 10.25];    


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
        mean(positions(13:15)) mean(positions(16:18)) mean(positions(19:21)) mean(positions(22:24))])
    set(gca,'xticklabel',{'1','2','3','4','5','6','7','8'})


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


    xlabel('Gradient set', 'fontsize', fontsize)
    ylabel('Relative error [%]', 'fontsize', fontsize);
    title([params_string(i,:)], 'fontsize', fontsize)
    
    %pause;
    w = 12;
    h = 6;
    set(gcf, 'PaperSize', [w h]);
    set(gcf, 'PaperPosition', [0 0 w h]);
    saveas(figure(i), fullfile(folder_pool, ['simexp3_' model '_' params_string(i,:)]), 'pdf');

end

%% MK/MD DKI, all sets, all tissues (groups of 3)

folder_pool = ['/Users/live/Documents/Masterthesis/MATLAB/DKI/Results_corrected/simulation_experiment_3/Pool/'];


fa_pool = [fa_DKI_WM fa_DKI_GM fa_DKI_CSF];
md_pool = [md_DKI_WM md_DKI_GM md_DKI_CSF];
rd_pool = [rd_DKI_WM rd_DKI_GM rd_DKI_CSF];
ad_pool = [ad_DKI_WM ad_DKI_GM ad_DKI_CSF];

mk_pool = [mk_DKI_WM(:,2:8) mk_DKI_GM(:,2:8) mk_DKI_CSF(:,2:8)];
rk_pool = [rk_DKI_WM(:,2:8) rk_DKI_GM(:,2:8) rk_DKI_CSF(:,2:8)];
ak_pool = [ak_DKI_WM(:,2:8) ak_DKI_GM(:,2:8) ak_DKI_CSF(:,2:8)];

model = 'DKI';
params_pool = [{ak_pool}]; 
params_string = ['AK'];

%group = [1 4 7 10 13 16 19 22 2 5 8 11 14 17 20 23 3 6 9 12 15 18 21 24];
%positions = [1 1.25 1.5 2.25 2.5 2.75 3.5 3.75 4 4.75 5 5.25 6 6.25 6.5 7.25 7.5 7.75 8.5 8.75 9 9.75 10 10.25];    

group = [1 4 7 10 13 16 19 22 2 5 8 11 14 17 20 23 3 6 9 12 15];
positions = [1 1.25 1.5 2.25 2.5 2.75 3.5 3.75 4 4.75 5 5.25 6 6.25 6.5 7.25 7.5 7.75 8.5 8.75 9];    

colorWM =  [0.98, 0.91, 0.71]; % white
colorGM = [0.41, 0.41, 0.41]; % gray
colorCSF = [0.34, 0.63, 0.83]; % blue

color = [{colorCSF},{colorGM},{colorWM},{colorCSF},{colorGM},{colorWM},{colorCSF},{colorGM},{colorWM},{colorCSF},{colorGM},{colorWM},{colorCSF},{colorGM},{colorWM},...
    {colorCSF},{colorGM},{colorWM},{colorCSF},{colorGM},{colorWM},{colorCSF},{colorGM},{colorWM}];


for i = 1:length(params_pool)
    x = params_pool{i};
    figure(i), boxplot(x, group, 'positions', positions);
    
    hold on
    plot(linspace(0,105), 0*ones(1,100), '-');
    hold off

    % Group name
    %set(gca,'xtick',[mean(positions(1:3)) mean(positions(4:6)) mean(positions(7:9)) mean(positions(10:12)) ...
    %    mean(positions(13:15)) mean(positions(16:18)) mean(positions(19:21)) mean(positions(22:24))]);
    %set(gca,'xticklabel',{'1','2','3','4','5','6','7','8',})    
    set(gca,'xtick',[mean(positions(1:3)) mean(positions(4:6)) mean(positions(7:9)) mean(positions(10:12)) ...
        mean(positions(13:15)) mean(positions(16:18)) mean(positions(19:21))]);
    set(gca,'xticklabel',{'2','3','4','5','6','7','8',})    


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


    xlabel('Gradient set', 'fontsize', fontsize)
    ylabel('Relative error [%]', 'fontsize', fontsize);
    title([params_string(i,:)], 'fontsize', fontsize)
    
    %pause;
    w = 12;
    h = 6;
    set(gcf, 'PaperSize', [w h]);
    set(gcf, 'PaperPosition', [0 0 w h]);
    saveas(figure(i), fullfile(folder_pool, ['simexp3_' model '_' params_string(i,:)]), 'pdf');

end