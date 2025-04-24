%% Modify as wanted
clear all
close all

strainFolder = 'CyaB'; % IMPORTANT!!! must match exactly the folder name, folder must contains "data" and "final_figure" folder, "data" must contain "Input_Profile_Analysis.xlsx"
fluo = 'mNeonGreen'; % IMPORTANT!!! must match exactly what was entered in BacStalk
% fluo = 'mScarletI'; % IMPORTANT!!! must match exactly what was entered in BacStalk

only_plot = 0;
load_name = 'Data_20241002_Strains_1756_2582_2588_2589_mNeonGreen'; % 'Data_Name' no .mat at the end
dir_data=strcat("/Volumes/Gani_sv_WS/bs_locprof_data_storage/",strainFolder);
dir_save_mat = strcat("/Volumes/Gani_sv_WS/bs_locprof_data_storage/",strainFolder,"/mat");
dir_save_graph = strcat("/Volumes/Gani_sv_WS/bs_locprof_results/",strainFolder);

% load functions
dir_func = '/Volumes/Gani_sv_WS/git/bs_locprof/Functions';
addpath(dir_func);

do_save = 1;% if 0 doesn't save graphs
do_violin = 1; % Note, doesn't work after using BacStalk. Restart Matlab!

%% Compute profiles and save
if ~only_plot
    tic
    clearvars -except strainFolder do_save only_plot do_violin fluo dir_save_mat dir_save_graph dir_data; close all; clc;

    % Field name of fluorescence profile data in the BacStalk Struct Dataset
    fluo_FieldNames=strcat('MedialAxisIntensity_',fluo);

    % Compute profile stats
    NbUpSamples=200; % Up sample the intensity profile to get similar length vectors to compare.
    k=linspace(0,1,NbUpSamples); %k is a ordered vector from 1 to 200
    lengthTH = 5;
    cutup = 1; % if 1 analyses cells BELOW "lengthTH"reshold, if 0 analyses cells above threshold
    CellLengthLiq=lengthTH; % in µm;
    CellLengthSol=lengthTH; % in µm;
    BootN=300;

    % load functions
    dir_func = '/Volumes/Gani_sv_WS/git/bs_locprof/Functions';
    addpath(dir_func);

    % Read Input file (Format of columns must be: Strain number - Strain name - date - condition) Note: condition names must all have the same number of characters
    [num,txt,~]=xlsread(strcat(dir_data, '/Input_Profile_Analysis.xlsx')); % must be located in 'data'
    Pil_numbers = num(:,1);
    Pil_numbers_unique = unique(Pil_numbers);
    dates = num(:,3); % read as a column vector
    dates_unique = unique(dates);
    Pil_types = txt(:,1); % read as a cell with one column
    Pil_types_unique = unique(txt(:,1));
    nbr_PTU = size(Pil_types_unique,1);
    conditions = txt(:,3); % read as a cell with one column
    conditions_unique = unique(conditions);
    nbr_conditions = size(conditions_unique,1);
    clear num txt

    bio_reps = cell(nbr_PTU,nbr_conditions);

    % extract info from BacStalk file and put into MegaStructure
    for type = 1:nbr_PTU
        Pil_type = Pil_types_unique{type};
        index_type = find(matches(Pil_types,Pil_type)==1);

        number_type = unique(Pil_numbers(index_type));
        dates_type = dates(index_type);

        conditions_type = conditions(index_type);
        conditions_type_unique = unique(conditions_type);
        nbr_conditions_type = size(conditions_type_unique,1);

        MegaStruct.(char(Pil_type)).full_name = strcat(num2str(number_type),'_',Pil_type);

        for cond = 1:nbr_conditions_type
           condition = conditions_type_unique(cond);
           index_condition = find(contains(conditions_type,condition)==1);
           dates_cond = dates_type(index_condition);
           index_conditions_unique = find(contains(conditions_unique, condition)==1);
           bio_reps{type,index_conditions_unique} = size(index_condition,1);

           dates_type_unique = unique(dates(index_type));
           MegaStruct.(char(Pil_type)).dates = dates_type_unique;
           Means = zeros(bio_reps{type,cond}, NbUpSamples);
           totalCells = 0;
           to_delete = zeros(1,bio_reps{type,cond});

           for rep = 1:bio_reps{type,index_conditions_unique}
               file_loc = strcat('/',num2str(dates_cond(rep)),'_DataSet_',num2str(number_type),'_',Pil_type,'_',condition);
               disp(strcat("File: ",file_loc{1}(2:end)));
               Path = strcat(dir_data, "/data/", file_loc,".mat");
                
               % Profiles                   
                   [Mean, Std, N, Profile, FluoMeans, CellWidth, CellLength, CellID, PixSize]=getMeanProfile(Path, k, NbUpSamples, fluo_FieldNames, CellLengthSol, cutup, fluo);
                   temp=trapz(Mean);
                   Mean=Mean/temp;

                   MegaStruct.(char(Pil_type)).(char(condition)).(char(strcat("BR_", num2str(dates_type_unique(rep))))).Mean = Mean;
                   MegaStruct.(char(Pil_type)).(char(condition)).(char(strcat("BR_", num2str(dates_type_unique(rep))))).Std = Std;
                   MegaStruct.(char(Pil_type)).(char(condition)).(char(strcat("BR_", num2str(dates_type_unique(rep))))).N = N;
                   MegaStruct.(char(Pil_type)).(char(condition)).(char(strcat("BR_", num2str(dates_type_unique(rep))))).Profile = Profile;
                   MegaStruct.(char(Pil_type)).(char(condition)).(char(strcat("BR_", num2str(dates_type_unique(rep))))).FluoMeans = FluoMeans;
                   MegaStruct.(char(Pil_type)).(char(condition)).(char(strcat("BR_", num2str(dates_type_unique(rep))))).CellWidth = CellWidth;
                   MegaStruct.(char(Pil_type)).(char(condition)).(char(strcat("BR_", num2str(dates_type_unique(rep))))).CellLength = CellLength;
                   MegaStruct.(char(Pil_type)).(char(condition)).(char(strcat("BR_", num2str(dates_type_unique(rep))))).CellID = CellID;

                   no_cells = isnan(mean(MegaStruct.(char(Pil_type)).(char(condition)).(char(strcat("BR_", num2str(dates_type_unique(rep))))).Mean));
                   if no_cells
                       to_delete(rep) = 1;
                   else
                       Means(rep,:) = MegaStruct.(char(Pil_type)).(char(condition)).(char(strcat("BR_", num2str(dates_type_unique(rep))))).Mean;
                   end

                   % Profiles condensed: Ratios
                   normalize = 1;
                   RefCytoNorm = 1; % Normalizes by taking the curve of cytoplasmic mNG expressed from a plasmid ad background signal
                   median_min = 1; % 1 = average middle fluorescence used as background, 0 = minimum middle fluorescence used as background
                   ratio_max = 0;
                   [PolarLoc, BootPolarLoc, MedianPolarLoc, CIPolarLoc, Polarisation, MedianPolarisation, NbCells, PoleInt]=polarRatios_MJK(Profile, CellWidth, CellLength, k, BootN, normalize, ratio_max, PixSize, median_min, RefCytoNorm);

                   MegaStruct.(char(Pil_type)).(char(condition)).(char(strcat("BR_", num2str(dates_type_unique(rep))))).PolarLoc = PolarLoc;
                   MegaStruct.(char(Pil_type)).(char(condition)).(char(strcat("BR_", num2str(dates_type_unique(rep))))).BootPolarLoc = BootPolarLoc;
                   MegaStruct.(char(Pil_type)).(char(condition)).(char(strcat("BR_", num2str(dates_type_unique(rep))))).MedianPolarLoc = MedianPolarLoc;
                   MegaStruct.(char(Pil_type)).(char(condition)).(char(strcat("BR_", num2str(dates_type_unique(rep))))).CIPolarLoc = CIPolarLoc;
                   MegaStruct.(char(Pil_type)).(char(condition)).(char(strcat("BR_", num2str(dates_type_unique(rep))))).Polarisation = Polarisation;
                   MegaStruct.(char(Pil_type)).(char(condition)).(char(strcat("BR_", num2str(dates_type_unique(rep))))).MedianPolarisation = MedianPolarisation;
                   MegaStruct.(char(Pil_type)).(char(condition)).(char(strcat("BR_", num2str(dates_type_unique(rep))))).NbCells = NbCells;
                   MegaStruct.(char(Pil_type)).(char(condition)).(char(strcat("BR_", num2str(dates_type_unique(rep))))).PoleInt = PoleInt;

                   totalCells = totalCells+NbCells;                 
           end
           Means(find(to_delete),:) = [];
           MegaStruct.(char(Pil_type)).(char(condition)).Means = Means;
           MegaStruct.(char(Pil_type)).(char(condition)).totalCells = totalCells;
        end    
    end

    % Save ALL
    clearvars -except fluo do_violin do_save only_plot dir_save_mat dir_save_graph bio_reps data_root_path dates conditions conditions_unique nbr_conditions dates_unique MegaStruct nbr_PTU NbUpSamples Pil_numbers Pil_numbers_unique Pil_types Pil_types_unique 
    if do_save
        save_name = 'Data';
        for segment = 1:size(dates_unique)
            save_name = strcat(save_name,'_',num2str(dates_unique(segment)));
        end
        save_name = strcat(save_name,'_Strains');
        for segment = 1:size(Pil_numbers_unique)
            save_name = strcat(save_name,'_',num2str(Pil_numbers_unique(segment)));
        end
        save(strcat(dir_save_mat,'/',save_name,'_',fluo,'.mat'));
    end
    toc
end

%% Plot Stuff
if only_plot
    close all
    do_save_onlyplot = do_save;
    clearvars -except fluo load_name strainFolder do_violin do_save_onlyplot dir_save_mat dir_save_graph
    load(strcat(dir_save_mat,'/',load_name,'.mat'));
    do_save = do_save_onlyplot;
    clear load_name
    tic
end

k=linspace(-1,1,NbUpSamples);
alpha=0.15;
linewidth=2;
dim=1;
y_limit=[0.003, 0.008];
filled_plot = 1;

colours = [236 32 36; 0 166 156; 139 197 63; 242 145 0; 29 112 183; 229 0 126];
if nbr_PTU > 6 | nbr_conditions > 6
    warning("WARNING! DANGER! There's not enough colours for the strains you want to plot, man!!!")
    return
end

nbr_columns = round((nbr_PTU+nbr_conditions)/2);

if matches(fluo, 'mNeonGreen')
    fluo = 'mNG';
elseif matches(fluo, 'mScarletI')
    fluo = 'mScI';
end

%% Plot Profiles

figure('units','normalized','outerposition',[0 0 1 1])
counter=0;

% plot for individual conditions
for cond = 1:nbr_conditions
    counter = counter+1;
    condition = conditions_unique{cond};
    subplot(2,nbr_columns,counter);
    leg = [];

    for type = 1:nbr_PTU
        Pil_type = Pil_types_unique{type};
                
        field_exists = isfield(MegaStruct.(char(Pil_type)),(char(condition)));
        if field_exists
            reps_cond = size(MegaStruct.(char(Pil_type)).(char(condition)).Means,1);
            total_cells = MegaStruct.(char(Pil_type)).(char(condition)).totalCells;
            leg = [leg,strcat(string(MegaStruct.(char(Pil_type)).full_name)," (reps: ",num2str(reps_cond),", cells: ",num2str(total_cells),")")];
            what2plot = MegaStruct.(char(Pil_type)).(char(condition)).Means;
            hold on
            plot(k, mean(what2plot,dim), 'Color', colours(type,:)/256,'LineWidth',linewidth);
        end
    end

    for type = 1:nbr_PTU
        Pil_type = Pil_types_unique{type};
        
        field_exists = isfield(MegaStruct.(char(Pil_type)),(char(condition)));
        if field_exists
            Means = MegaStruct.(char(Pil_type)).(char(condition)).Means;
            top = mean(Means,dim)+ std(Means);
            bottom = mean(Means,dim)- std(Means);
            k2 = [k, fliplr(k)];
            if filled_plot
                inBetween = [top, fliplr(bottom)];
                fill(k2, inBetween, colours(type,:)/256,'FaceAlpha',alpha, 'EdgeColor', 'none');
            else
                plot(k, top, k, bottom, 'Color', colours(type,:)/256, 'LineStyle', '--','LineWidth',linewidth-1);
            end  
        end
    end

    legend(leg, 'Location', 'Best', 'Interpreter', 'none');
    ylim(y_limit);
    xlabel('Normalized distance from midcell');
    ylabel(strcat('Normalized fluorescence',' (',fluo,')'));
    title(condition, 'Interpreter', 'none');   

end

% plot for individual strains
load('ReferenceCytoplasmicCurve.mat');
for type = 1:nbr_PTU
    counter = counter+1;
    subplot(2,nbr_columns,counter);
    Pil_type = Pil_types_unique{type};
    leg = []; 

    for cond = 1:nbr_conditions
        condition = conditions_unique{cond};
        
        field_exists = isfield(MegaStruct.(char(Pil_type)),(char(condition)));
        if field_exists
            reps_cond = size(MegaStruct.(char(Pil_type)).(char(condition)).Means,1);
            total_cells = MegaStruct.(char(Pil_type)).(char(condition)).totalCells;
            leg = [leg;strcat(string(condition)," (reps: ",num2str(reps_cond),", cells: ",num2str(total_cells),")")];        
            what2plot = MegaStruct.(char(Pil_type)).(char(condition)).Means;
            hold on
            plot(k, mean(what2plot,dim), 'Color', colours(cond,:)/256,'LineWidth',linewidth); 
        end
    end
    for cond = 1:nbr_conditions
        condition = conditions_unique{cond};
        
        field_exists = isfield(MegaStruct.(char(Pil_type)),(char(condition)));
        if field_exists
            Means = MegaStruct.(char(Pil_type)).(char(condition)).Means;
            top = mean(Means,dim)+ std(Means);
            bottom = mean(Means,dim)- std(Means);
            k2 = [k, fliplr(k)];
                if filled_plot
                    inBetween = [top, fliplr(bottom)];
                    fill(k2, inBetween, colours(cond,:)/256,'FaceAlpha',alpha, 'EdgeColor', 'none');
                else
                    plot(k, top, k, bottom, 'Color', colours(cond,:)/256, 'LineStyle', '--','LineWidth',linewidth-1);
                end  
        end
    end
    plot(k, ReferenceCytoplasmicCurve, '--' ,'Color', [128 128 128]/256,'LineWidth',1);
    legend(leg, 'Location', 'Best', 'Interpreter', 'none');
    ylim(y_limit);
    xlabel('Normalized distance from midcell');
    ylabel(strcat('Normalized fluorescence',' (',fluo,')'));
    title(string(MegaStruct.(char(Pil_type)).full_name), 'Interpreter', 'none'); 
end

% save
if do_save
    save_name = '';
    for segment = 1:size(dates_unique)
        save_name = strcat(save_name,num2str(dates_unique(segment)),'_');
    end
    save_name = strcat(save_name,'Strains_');
    for segment = 1:size(Pil_numbers_unique)
        save_name = strcat(save_name,num2str(Pil_numbers_unique(segment)),'_');
    end
    save_name = strcat(save_name,'Profiles','_',fluo);
    saveas(gcf,strcat(dir_save_graph,'/fig_files/',save_name,'.fig'));
    saveas(gcf,strcat(dir_save_graph,'/svg_files/',save_name,'.svg'));
    saveas(gcf,strcat(dir_save_graph,'/',save_name,'.jpg'));  
end
%% Plot total fluorescence

figure('units','normalized','outerposition',[0 0 1 1])
counter=0;

y_max = 300;

% plot for individual conditions
for cond = 1:nbr_conditions
    counter = counter+1;
    condition = conditions_unique{cond};
    subplot(2,nbr_columns,counter);
    xlabels = [];   
    position = 0;

    for type = 1:nbr_PTU
        Pil_type = Pil_types_unique{type};
        xlabels = [xlabels;string(MegaStruct.(char(Pil_type)).full_name)];
        what2plot = [];
        position = position+1;
        
        field_exists = isfield(MegaStruct.(char(Pil_type)),(char(condition)));
        if field_exists
            for rep = 1:bio_reps{type,cond}
                BR_name = strcat("BR_", num2str(MegaStruct.(char(Pil_type)).dates(rep)));
                what2plot = [what2plot, mean(MegaStruct.(char(Pil_type)).(char(condition)).(char(BR_name)).FluoMeans)];   
            end
            hold on
            what2plot = rmmissing(what2plot);
            plot(position, what2plot, "o", "MarkerSize",6.5, "Linewidth",1.5, 'Color', colours(type,:)/256)
            plot([position-0.1 position+0.1], [mean(what2plot) mean(what2plot)],"-k","MarkerSize",6.5,"Linewidth",1.5) 
        end

    end
    axis([0 position+1 0 y_max])
    set(gca, 'XTickLabel',["";xlabels;""], 'Fontsize',10, 'TickLabelInterpreter','none')
    xticks([0:position+1])
    xtickangle(15)
    ylabel(strcat('Mean fluorescence',' (',fluo,')'));
    title(condition, 'Interpreter', 'none');   
end


% plot for individual strains
for type = 1:nbr_PTU
    counter = counter+1;
    subplot(2,nbr_columns,counter);
    Pil_type = Pil_types_unique{type};
    position = 0;
    xlabels = [];

    for cond = 1:nbr_conditions
        condition = conditions_unique{cond};
        xlabels = [xlabels; condition];
        what2plot = [];
        position = position+1;
        
        field_exists = isfield(MegaStruct.(char(Pil_type)),(char(condition)));
        if field_exists            
            for rep = 1:bio_reps{type,cond}
                BR_name = strcat("BR_", num2str(MegaStruct.(char(Pil_type)).dates(rep)));
                what2plot = [what2plot, mean(MegaStruct.(char(Pil_type)).(char(condition)).(char(BR_name)).FluoMeans)];   
            end
            hold on
            what2plot = rmmissing(what2plot);
            plot(position, what2plot, "o", "MarkerSize",6.5, "Linewidth",1.5, 'Color', colours(type,:)/256)
            plot([position-0.1 position+0.1], [mean(what2plot) mean(what2plot)],"-k","MarkerSize",6.5,"Linewidth",1.5)
        end
    end
    axis([0 position+1 0 y_max])
    set(gca, 'XTickLabel',["";xlabels;""], 'Fontsize',10, 'TickLabelInterpreter','none')
    xticks([0:position+1])
    xtickangle(15)
    ylabel(strcat('Mean fluorescence',' (',fluo,')'));
    title(string(MegaStruct.(char(Pil_type)).full_name), 'Interpreter', 'none');   
end

% save
if do_save
    save_name = '';
    for segment = 1:size(dates_unique)
        save_name = strcat(save_name,num2str(dates_unique(segment)),'_');
    end
    save_name = strcat(save_name,'Strains_');
    for segment = 1:size(Pil_numbers_unique)
        save_name = strcat(save_name,num2str(Pil_numbers_unique(segment)),'_');
    end
    save_name = strcat(save_name,'MeanFluo','_',fluo);
    saveas(gcf,strcat(dir_save_graph,'/fig_files/',save_name,'.fig'));
    saveas(gcf,strcat(dir_save_graph,'/svg_files/',save_name,'.svg'));
    saveas(gcf,strcat(dir_save_graph,'/',save_name,'.jpg'));
end
%% Plot Polar Locaization Index

figure('units','normalized','outerposition',[0 0 1 1])
counter=0;

y_min = 0;
y_max = 0.6;
vioscale = 0.15;

% plot for individual conditions
for cond = 1:nbr_conditions
    counter = counter+1;
    condition = conditions_unique{cond};
    subplot(2,nbr_columns,counter);
    xlabels = [];   
    position = 0;

    for type = 1:nbr_PTU
        Pil_type = Pil_types_unique{type};
        xlabels = [xlabels;string(MegaStruct.(char(Pil_type)).full_name)];
        what2plot = [];
        vioplot = [];
        position = position+1;
        
        field_exists = isfield(MegaStruct.(char(Pil_type)),(char(condition)));
        if field_exists
            for rep = 1:bio_reps{type,cond}
                BR_name = strcat("BR_", num2str(MegaStruct.(char(Pil_type)).dates(rep)));
                what2plot = [what2plot, MegaStruct.(char(Pil_type)).(char(condition)).(char(BR_name)).MedianPolarLoc];
                vioplot = [vioplot, MegaStruct.(char(Pil_type)).(char(condition)).(char(BR_name)).PolarLoc];                
            end
            hold on
            what2plot = rmmissing(what2plot);
            if do_violin
                violin(position,vioplot,'facealpha',0.42,'linewidth',0.5,'facecolor',colours(type,:)/256,'style',2,'side','right','scaling',vioscale)
            end
            plot(position, what2plot, "o", "MarkerSize",6.5, "Linewidth",1.5, 'Color', colours(type,:)/256)
            plot([position-0.1 position+0.1], [mean(what2plot) mean(what2plot)],"-k","MarkerSize",6.5,"Linewidth",1.5) 
        end

    end
    axis([0 position+1 y_min y_max])
    set(gca, 'XTickLabel',["";xlabels;""], 'Fontsize',10, 'TickLabelInterpreter','none')
    xticks([0:position+1])
    xtickangle(15)
    ylabel(strcat('Median Polar Localization Index',' (',fluo,')'));
    title(condition, 'Interpreter', 'none');   
end

% plot for individual strains
for type = 1:nbr_PTU
    counter = counter+1;
    subplot(2,nbr_columns,counter);
    Pil_type = Pil_types_unique{type};
    position = 0;
    xlabels = [];

    for cond = 1:nbr_conditions
        condition = conditions_unique{cond};
        xlabels = [xlabels; condition];
        what2plot = [];
        vioplot = [];
        position = position+1;
        
        field_exists = isfield(MegaStruct.(char(Pil_type)),(char(condition)));
        if field_exists            
            for rep = 1:bio_reps{type,cond}
                BR_name = strcat("BR_", num2str(MegaStruct.(char(Pil_type)).dates(rep)));
                what2plot = [what2plot, MegaStruct.(char(Pil_type)).(char(condition)).(char(BR_name)).MedianPolarLoc];
                vioplot = [vioplot, MegaStruct.(char(Pil_type)).(char(condition)).(char(BR_name)).PolarLoc];
            end
            hold on
            what2plot = rmmissing(what2plot);
            if do_violin
                violin(position,vioplot,'facealpha',0.42,'linewidth',0.5,'facecolor',colours(type,:)/256,'style',2,'side','right','scaling',vioscale)
            end
            plot(position, what2plot, "o", "MarkerSize",6.5, "Linewidth",1.5, 'Color', colours(type,:)/256)
            plot([position-0.1 position+0.1], [mean(what2plot) mean(what2plot)],"-k","MarkerSize",6.5,"Linewidth",1.5)
        end
    end
    axis([0 position+1 y_min y_max])
    set(gca, 'XTickLabel',["";xlabels;""], 'Fontsize',10, 'TickLabelInterpreter','none')
    xticks([0:position+1])
    xtickangle(15)
    ylabel(strcat('Median Polar Localization Index',' (',fluo,')'));
    title(string(MegaStruct.(char(Pil_type)).full_name), 'Interpreter', 'none');   
end

% save
if do_save
    save_name = '';
    for segment = 1:size(dates_unique)
        save_name = strcat(save_name,num2str(dates_unique(segment)),'_');
    end
    save_name = strcat(save_name,'Strains_');
    for segment = 1:size(Pil_numbers_unique)
        save_name = strcat(save_name,num2str(Pil_numbers_unique(segment)),'_');
    end
    save_name = strcat(save_name,'PolarLocalizationIndex','_',fluo);
    saveas(gcf,strcat(dir_save_graph,'/fig_files/',save_name,'.fig'));
    saveas(gcf,strcat(dir_save_graph,'/svg_files/',save_name,'.svg'));
    saveas(gcf,strcat(dir_save_graph,'/',save_name,'.jpg'));
end

%% Plot Symmetry Index (Polarization)

figure('units','normalized','outerposition',[0 0 1 1])
counter=0;

y_min = 0.5;
y_max = 0.6;

% plot for individual conditions
for cond = 1:nbr_conditions
    counter = counter+1;
    condition = conditions_unique{cond};
    subplot(2,nbr_columns,counter);
    xlabels = [];   
    position = 0;

    for type = 1:nbr_PTU
        Pil_type = Pil_types_unique{type};
        xlabels = [xlabels;string(MegaStruct.(char(Pil_type)).full_name)];
        what2plot = [];
        position = position+1;
        
        field_exists = isfield(MegaStruct.(char(Pil_type)),(char(condition)));
        if field_exists
            for rep = 1:bio_reps{type,cond}
                BR_name = strcat("BR_", num2str(MegaStruct.(char(Pil_type)).dates(rep)));
                what2plot = [what2plot, MegaStruct.(char(Pil_type)).(char(condition)).(char(BR_name)).MedianPolarisation];   
            end
            hold on
            what2plot = rmmissing(what2plot);
            plot(position, what2plot, "o", "MarkerSize",6.5, "Linewidth",1.5, 'Color', colours(type,:)/256)
            plot([position-0.1 position+0.1], [mean(what2plot) mean(what2plot)],"-k","MarkerSize",6.5,"Linewidth",1.5) 
        end

    end
    axis([0 position+1 y_min y_max])
    set(gca, 'XTickLabel',["";xlabels;""], 'Fontsize',10, 'TickLabelInterpreter','none')
    xticks([0:position+1])
    xtickangle(15)
    ylabel(strcat('Median Symmetry Index',' (',fluo,')'));
    title(condition, 'Interpreter', 'none');   
end

% plot for individual strains
for type = 1:nbr_PTU
    counter = counter+1;
    subplot(2,nbr_columns,counter);
    Pil_type = Pil_types_unique{type};
    position = 0;
    xlabels = [];

    for cond = 1:nbr_conditions
        condition = conditions_unique{cond};
        xlabels = [xlabels; condition];
        what2plot = [];
        position = position+1;
        
        field_exists = isfield(MegaStruct.(char(Pil_type)),(char(condition)));
        if field_exists            
            for rep = 1:bio_reps{type,cond}
                BR_name = strcat("BR_", num2str(MegaStruct.(char(Pil_type)).dates(rep)));
                what2plot = [what2plot, MegaStruct.(char(Pil_type)).(char(condition)).(char(BR_name)).MedianPolarisation];   
            end
            hold on
            what2plot = rmmissing(what2plot);
            plot(position, what2plot, "o", "MarkerSize",6.5, "Linewidth",1.5, 'Color', colours(type,:)/256)
            plot([position-0.1 position+0.1], [mean(what2plot) mean(what2plot)],"-k","MarkerSize",6.5,"Linewidth",1.5)
        end
    end
    axis([0 position+1 y_min y_max])
    set(gca, 'XTickLabel',["";xlabels;""], 'Fontsize',10, 'TickLabelInterpreter','none')
    xticks([0:position+1])
    xtickangle(15)
    ylabel(strcat('Median Symmetry Index',' (',fluo,')'));
    title(string(MegaStruct.(char(Pil_type)).full_name), 'Interpreter', 'none');   
end

% save
if do_save
    save_name = '';
    for segment = 1:size(dates_unique)
        save_name = strcat(save_name,num2str(dates_unique(segment)),'_');
    end
    save_name = strcat(save_name,'Strains_');
    for segment = 1:size(Pil_numbers_unique)
        save_name = strcat(save_name,num2str(Pil_numbers_unique(segment)),'_');
    end
    save_name = strcat(save_name,'SymmetryIndex','_',fluo);
    saveas(gcf,strcat(dir_save_graph,'/fig_files/',save_name,'.fig'));
    saveas(gcf,strcat(dir_save_graph,'/svg_files/',save_name,'.svg'));
    saveas(gcf,strcat(dir_save_graph,'/',save_name,'.jpg'));
end

%% Plot Polar Intensities against Cell Length

% figure('units','normalized','outerposition',[0 0 1 1])
% counter=0;
% 
% x_max = 10;
% y_min = 0;
% y_max = 0.016;
% 
% % plot for individual strains - Max int
% for type = 1:nbr_PTU
%     counter = counter+1;
%     subplot(2,nbr_columns,counter);
%     Pil_type = Pil_types_unique{type};
%     
%     for cond = 1:nbr_conditions
%         condition = conditions_unique{cond};
%         what2plot = [];
%         x_vals = [];
%         
%         field_exists = isfield(MegaStruct.(char(Pil_type)),(char(condition)));
%         if field_exists            
%             for rep = 1:bio_reps{type,cond}
%                 BR_name = strcat("BR_", num2str(MegaStruct.(char(Pil_type)).dates(rep)));
%                 x_vals = [x_vals, MegaStruct.(char(Pil_type)).(char(condition)).(char(BR_name)).PoleInt(1,:)];  
%                 what2plot = [what2plot, MegaStruct.(char(Pil_type)).(char(condition)).(char(BR_name)).PoleInt(2,:)];   
%             end
%             hold on
%             what2plot = rmmissing(what2plot);
%             plot(x_vals, what2plot, ".", "MarkerSize",6.5, "Linewidth",1.5, 'Color', colours(type,:)/256)
%             stdev = std(what2plot);
%             text(1,0.005,num2str(stdev))
%          end
%     end
%     axis([0 x_max y_min y_max])
%     xlabel('Cell Length (µm)');
%     ylabel('Absolute Polar Intensity');
%     title(strcat(string(MegaStruct.(char(Pil_type)).full_name)," - Bright Pole"), 'Interpreter', 'none');   
% end
% 
% % plot for individual strains - Min int
% for type = 1:nbr_PTU
%     counter = counter+1;
%     subplot(2,nbr_columns,counter);
%     Pil_type = Pil_types_unique{type};
%     
%     for cond = 1:nbr_conditions
%         condition = conditions_unique{cond};
%         what2plot = [];
%         x_vals = [];
%         
%         field_exists = isfield(MegaStruct.(char(Pil_type)),(char(condition)));
%         if field_exists            
%             for rep = 1:bio_reps{type,cond}
%                 BR_name = strcat("BR_", num2str(MegaStruct.(char(Pil_type)).dates(rep)));
%                 x_vals = [x_vals, MegaStruct.(char(Pil_type)).(char(condition)).(char(BR_name)).PoleInt(1,:)];  
%                 what2plot = [what2plot, MegaStruct.(char(Pil_type)).(char(condition)).(char(BR_name)).PoleInt(3,:)];   
%             end
%             hold on
%             what2plot = rmmissing(what2plot);
%             plot(x_vals, what2plot, ".", "MarkerSize",6.5, "Linewidth",1.5, 'Color', colours(type,:)/256)
%             stdev = std(what2plot);
%             text(1,0.005,num2str(stdev))
%          end
%     end
%     axis([0 x_max y_min y_max])
%     xlabel('Cell Length (µm)');
%     ylabel('Absolute Polar Intensity');
%     title(strcat(string(MegaStruct.(char(Pil_type)).full_name)," - Dim Pole"), 'Interpreter', 'none');  
% end
% 
% % save
% if do_save
%     save_name = '';
%     for segment = 1:size(dates_unique)
%         save_name = strcat(save_name,num2str(dates_unique(segment)),'_');
%     end
%     save_name = strcat(save_name,'Strains_');
%     for segment = 1:size(Pil_numbers_unique)
%         save_name = strcat(save_name,num2str(Pil_numbers_unique(segment)),'_');
%     end
%     save_name = strcat(save_name,'PoleIntensities');
%     % saveas(gcf,strcat(dir_save_graph,'/fig_files/',save_name,'.fig'));
%     % saveas(gcf,strcat(dir_save_graph,'/svg_files/',save_name,'.svg'));
%     saveas(gcf,strcat(dir_save_graph,'/',save_name,'.jpg'));
% end
toc