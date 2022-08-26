% Open corresponding analysis file Data_???_StrainNumber.mat and copy strain name and corresponding dates from MegaStruct
% Copy from export to us for SuperPlotsofData https://huygens.science.uva.nl/SuperPlotsOfData/
% Paste and select tab separated
strainName = 'mNG_PilH_D52A_cyaB';
dates = [20211019;20211020];
% conditions = {'sol0h'; 'sol1h'; 'sol2h'};
NumbrDates = size(dates,1);
NumbrConditions = size(conditions_unique,1);
% data_0h = cell(NumbrDates,2);
% toDel0h = zeros(NumbrDates,1);
% data_1h = cell(NumbrDates,2);
% toDel1h = zeros(NumbrDates,1);
% data_2h = cell(NumbrDates,2);
% toDel2h = zeros(NumbrDates,1);
header = cell(1,3);
header{1,1} = 'time on surface';
header{1,2} = 'Polar Localization Index';
header{1,3} = 'Asymmetry Index';
header{1,4} = 'Date';
header{1,5} = 'Strain';
export = [];
for c=1:NumbrConditions
    data = cell(NumbrDates,4);
    noSample = zeros(NumbrDates,1);
    condition = conditions_unique{c};  
    for i=1:NumbrDates
        name = strcat("BR_",num2str(dates(i)));
        fieldExists = isfield(MegaStruct.(char(strainName)).(char(condition)),(char(name)));
        if fieldExists    
            data{i,1} = condition;
            data{i,2} = MegaStruct.(char(strainName)).(char(condition)).(char(name)).MeanRat;
            data{i,3} = MegaStruct.(char(strainName)).(char(condition)).(char(name)).MeanPolarisation;
            data{i,4} = dates(i);
            data{i,5} = strainName;
        else
            noSample(i,1) = 1;
        end
    end
    toDel = find(noSample);
    data(toDel,:) = [];
    export = [export; data];
end
export = [header; export];
save_dir = 'H:\Marco_WS\Localization Profiles Automated Index\export for statistical tests\';
writecell(export,strcat(save_dir,strainName,'.csv'));