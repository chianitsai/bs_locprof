load('H:\Marco_WS\Localization Profiles Automated Index\Latest script versions\functions\ReferenceCytoplasmicCurve.mat');
% Theoretical profile with 100% signal at the poles
profile = zeros(1,200);
% RangePoles = 0.15; % which percent of the normalized cell length is considered the pole area; normally calcluated based on the width and length of cells
RangePoles = [0.05, 0.1, 0.15, 0.2, 0.25, 0.3];
PolarRatio = zeros(2,size(RangePoles,2));
for PZ=1:size(RangePoles,2)
   
    StartRange = 200*RangePoles(PZ);
    EndRange = 200*(1-RangePoles(PZ))+1;
    for i=1:200
        if i<=StartRange | i>=EndRange
            profile(i) = 1/(200*RangePoles(PZ)*2);
        else
            profile(i) = 0;
        end
    end

    Midline = ReferenceCytoplasmicCurve;
    AreaDiffusedStart=trapz(Midline(1:StartRange)); % calculates the theoretical mean fluorescence at pole 1 if signal was diffused everywhere
    AreaDiffusedEnd=trapz(Midline(EndRange:end)); % calculates the theoretical mean fluorescence at pole 2 if signal was diffused everywhere

    AreaStart=trapz(profile(1:StartRange)); % integrated fluorescence pole 1
    AreaEnd=trapz(profile(EndRange:end)); % integrated fluorescence pole 2

    I1=(AreaStart-AreaDiffusedStart); % corrected fluorescence pole 1 (substract theoretical diffused fluorescence)
    I2=(AreaEnd-AreaDiffusedEnd); % corrected fluorescence pole 2 (substract theoretical diffused fluorescence)
    
    PolarRatio(1,PZ) = RangePoles(PZ)*100;
    PolarRatio(2,PZ) = ((I1+I2)/(AreaStart+AreaEnd)); % in other words: PolarRatio = corrected polar fluorescence devided by total polar fluorescence

end       