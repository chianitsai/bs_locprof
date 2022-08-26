function [ PolarLoc, BootPolarLoc, MedianPolarLoc, CIPolarLoc, Polarisation, MedianPolarisation, NbCells, PoleInt] = polarRatios2( Profiles, CellWidths, CellLengths, k, BootN, normalize, ratio_max, PixSize, mean_min, RefCytoNorm)

    load('functions\ReferenceCytoplasmicCurve.mat');

    NbCells=size(Profiles,1);
    todel=[];
    for cell=1:NbCells
        if CellWidths(cell)>=CellLengths(cell)
            todel = [todel,cell];
        end
    end
    disp(strcat("Perfectly round or just weird cells: ", num2str(size(todel,2))));
    Profiles(todel,:)=[];
    CellWidths(todel)=[];
    CellLengths(todel)=[];

    n=0;
    NbCells=size(Profiles,1);
    dim=floor(sqrt(NbCells))+1;

    PolarRatio=zeros(1,NbCells);
    Polarism=zeros(1,NbCells);
    PoleInt=zeros(3,NbCells);
    for cell=1:NbCells
        %%%
        RangePoles=(CellWidths(cell))/(CellLengths(cell))*0.5; % calculates what % of the cell length is defined as the pole (cell width / cell length *0.5)

        StartRange=min(find(k>=(RangePoles))); % defines endpoint of pole region 1 in k linspace vector
        EndRange=max(find(k<=(1-RangePoles))); % defines endpoint of pole region 2 in k linspace vector

        profile=Profiles(cell,:);
        
        if RefCytoNorm
            Midline = ReferenceCytoplasmicCurve;
        else
            if mean_min
                MiddleFluorescence=median(profile(StartRange:EndRange)); % mean fluorescence of the non-polar region 
            else
                MiddleFluorescence=min(profile(StartRange:EndRange)); % min fluorescence of the non-polar region
            end
            SDMiddleFluorescence=std(profile(StartRange:EndRange)); % std of fluorescence of the non-polar region
            Midline=ones(1,length(k))*MiddleFluorescence; % defines a linspace like k but with mean fluorescence of the non-polar region as values
        end
        
        AreaDiffusedStart=trapz(Midline(1:StartRange)); % calculates the theoretical mean fluorescence at pole 1 if signal was diffused everywhere
        AreaDiffusedEnd=trapz(Midline(EndRange:end)); % calculates the theoretical mean fluorescence at pole 2 if signal was diffused everywhere

        AreaStart=trapz(profile(1:StartRange)); % integrated fluorescence pole 1
        AreaEnd=trapz(profile(EndRange:end)); % integrated fluorescence pole 2
        MaxStart=maxN(profile(1:StartRange),1);
        MaxEnd=maxN(profile(EndRange:end),1);
        
        % Polar Localization Index
        if normalize % normalizes the fluorescence at the poles because of lower polar signal artefact: If polar fluorescence below theoretical fluorescence estimated from cytoplasmic fluorescence (Midline) it sets it to 0 instead of below 0.
            if (((AreaStart-AreaDiffusedStart)>0)) % && (MaxStart>(MiddleFluorescence+2*SDMiddleFluorescence))) % first part checks if polar fluorescence > cytoplasmic fluorescence. not sure what second part does.
                I1=(AreaStart-AreaDiffusedStart); % corrected fluorescence pole 1 (substract theoretical diffused fluorescence)
            else
                I1=0;
            end
            if (((AreaEnd-AreaDiffusedEnd)>0)) % && (MaxEnd>(MiddleFluorescence+2*SDMiddleFluorescence)))
                I2=(AreaEnd-AreaDiffusedEnd); % corrected fluorescence pole 2 (substract theoretical diffused fluorescence)
            else
                I2=0;
            end
        else
            I1=(AreaStart-AreaDiffusedStart); % corrected fluorescence pole 1 (substract theoretical diffused fluorescence)
            I2=(AreaEnd-AreaDiffusedEnd); % corrected fluorescence pole 2 (substract theoretical diffused fluorescence)
        end
        
        PolarRatio(cell) = ((I1+I2)/(AreaStart+AreaEnd)); % in other words: PolarRatio = corrected polar fluorescence devided by total polar fluorescence
        
        % Symmetry Index
        if ratio_max
            I1=max(profile(1:StartRange));
            I2=max(profile(EndRange:end));
        else
            I1=(AreaStart);
            I2=(AreaEnd);
        end

        if(I1>=I2)
            Imax=I1;
            Imin=I2;
        else
            Imax=I2;
            Imin=I1;
        end
        if((I1+I2)>0)
            Polarism(cell) = Imax/(I1+I2); % in other words; Plarism = intensity at brighter pole devided by intensity of both poles
        else
            Polarism(cell) = 0.5; % if no polar fluorescence at all, set it to 0.5 (meaning same on both poles)
            n=n+1;
        end
        if (I2<0) % not sure what this does
            disp(strcat('I1 :', num2str(I1), ', I2 :', num2str(I2), ', I1+I2 :', num2str(I1+I2)));
            disp(strcat('PolarRatio :', num2str(PolarRatio(cell))));
            disp(strcat('MaxStart :', num2str(MaxStart), ' MaxEnd :', num2str(MaxEnd), ' Middle fluorescence :', num2str(MiddleFluorescence), ' Middle fluorescence SD:', num2str(SDMiddleFluorescence)));
        end
        % save cell length and Intensities of dim and bright pole in array
        PoleInt(1,cell) = CellLengths(cell)*PixSize;
        PoleInt(2,cell) = Imax;
        PoleInt(3,cell) = Imin;

    end
    NbCellsDiffused=n;
    disp(strcat('Total cells :', num2str(NbCells)));
    PolarLoc = PolarRatio(~isnan(PolarRatio));
    Polarisation = Polarism(~isnan(Polarism));
    if ~isempty(PolarRatio) | ~isempty(Polarism)
        BootPolarLoc =  bootstrp(BootN, @(x) median(x), PolarRatio);    
        MedianPolarLoc = median(PolarLoc);
        MedianPolarisation = median(Polarisation);
        CIPolarLoc(1,:) = ci_percentile(PolarLoc);
    else
        BootPolarLoc = [];    
        MedianPolarLoc = [];  
        MedianPolarisation = [];  
        CIPolarLoc = [];  
    end
    disp(" ");
end

