%% Identifying the type of gaits for a given solution or solution branch, return the name of gait and the color for visualization
function [gait,abbr, color_plot, linetype] = Gait_Identification_Asym(results)
    if size(results,2)<15
        [gait,abbr, color_plot, linetype] = Type_of_Gait(results(:,1));
    else
        downsample_rate = fix(size(results,2)/15);
        indices = downsample(1:size(results,2),downsample_rate);
        count = 1;
        
        [gait_last,abbr_current, color_current, linetype_current] = Type_of_Gait(results(:,indices(1)));

        for i = 2:length(indices)
            X = results(:,indices(i));
            [gait_current,abbr_current, color_current, linetype_current] = Type_of_Gait(X);

            if string(gait_current)==string(gait_last)
                count = count +1;
                if count>6
                    gait = gait_current;
                    abbr = abbr_current;
                    color_plot = color_current;
                    linetype = linetype_current;
                    return
                end
            end

            if i == length(indices)
                    disp('Gait identification error for current branch.')
                    disp('Using gait of current solution for representative.')
                    [gait,abbr, color_plot, linetype] = Type_of_Gait(X);
            end

            gait_last = gait_current;
        end

    end
end

%% Gait identifying function
function [gait,abbr, color_plot, linetype] = Type_of_Gait(X)
    
    gait = [];
    abbr = [];
    color_plot = [];
    linetype = [];


    threshold = 1e-6; % threshold for identifying 'equal'

    if abs(X(14)-X(18))<threshold && abs(X(16)-X(20))<threshold && abs(X(14)-X(16))<threshold && abs(X(16)-X(20))<threshold
        type1 = 'Pronking';
        abbr1 = 'PF';
        color_plot = [0 0.4470 0.7410];
    elseif abs(X(14)-X(18))<threshold && abs(X(16)-X(20))<threshold
        type1 = 'Bounding';
        abbr1 = 'B';
        color_plot = [0.8500 0.3270 0.0980];
    elseif abs(X(14)-X(18))<threshold 
        type1 = 'Half-Bounding with Front Legs Spread';
        abbr1 = 'F';
        color_plot = [0.4660 0.6740 0.1880];
    elseif abs(X(16)-X(20))<threshold 
        type1 = 'Half-Bounding with Hind Legs Spread';
        abbr1 = 'H';
        color_plot = [0.9290 0.6940 0.1270];
    else
        type1 = 'Galloping';
        abbr1 = 'G';
        color_plot = [0.4940 0.1840 0.5560];
    end
    
    nof = FlightPhaseCheck(X);
    MidstanceDiff = CalMidstance(X);
    % if ~(nof==1) && ( MidstanceDiff>0.3 && MidstanceDiff<0.70)
    if  MidstanceDiff>0.3 && MidstanceDiff<0.70
        type2 = '_with Additional Flight Phases';
        abbr2 = '2';
        linetype = '--';
    elseif abs(X(5))<1e-9
        type2 = '';
        abbr2 = '';
        linetype = '-';
    elseif X(5)>0
        type2 = '_with Gathered Suspension';
        abbr2 = 'G';
        linetype = '-';
    elseif  X(5)<0
        type2 = '_with Extended Suspension';
        abbr2 = 'E';
        linetype = ':';
    end
    
    gait = string(type1) + string(type2);
    abbr = string(abbr1) + string(abbr2);

end

%% Flight phase check
function nof = FlightPhaseCheck(X)
   X(14:22) = round(X(14:22),4);
   % draw time line
   dt = 1e-4;
  
   % extract the stance time of each leg
   if X(14)<X(15)
      bls = X(14):dt:X(15);
   else
      bls = [0:dt:X(15)  X(14):dt:X(22)];
   end 
   if X(16)<X(17)
      fls = X(16):dt:X(17);
   else
      fls = [0:dt:X(17)  X(16):dt:X(22)];
   end  
   if X(18)<X(19)
      brs = X(18):dt:X(19);
   else
      brs = [0:dt:X(19)  X(18):dt:X(22)];
   end 
   if X(20)<X(21)
      frs = X(20):dt:X(21);
   else
      frs = [0:dt:X(21)  X(20):dt:X(22)];
   end
   % extract total stance time
   stance = unique(sort(round([bls fls brs frs],4)));
   
   diffstance = round(diff(stance)/dt);
   nof = length(find(~(diffstance==1)));
   
   if X(14)<X(15) && X(16)<X(17) && X(18)<X(19) && X(20)<X(21)
       nof = nof +1;
   end
   
    % for bounding gait, if the flight phase is shorter than 1/5 of the
    % stride, the number of flight phase is considered as 1
    if abs(X(14)-X(18))<1e-6 && abs(X(16)-X(20))<1e-6 
        if sign(X(5))==1 && X(16)-X(15)<0.2*X(22)
            nof = 1;
        end
        if sign(X(5))==-1 && X(14)-X(17)<0.2*X(22)
            nof = 1;
        end
    end
end

%% Midstance Diff Check
function MidstanceDiff = CalMidstance(X)
    % calculate mid stance of left hind
    if X(14)<X(15)
        lh_ms = (X(14)+X(15))/2/X(22);
    elseif X(14)>X(15)
        lh_ms = (X(14)+X(15)+X(22))/2/X(22);
    end
    
    % calculate mid stance of right hind
    if X(18)<X(19)
        rh_ms = (X(18)+X(19))/2/X(22);
    elseif X(18)>X(19)
        rh_ms = (X(18)+X(19)+X(22))/2/X(22);
    end

    % calculate midstance of hind leg pair
    if norm(lh_ms-rh_ms)<0.5
        hms = (lh_ms+rh_ms)/2;
    else
        hms = (lh_ms + rh_ms + 1)/2;
    end
        
    if hms>1
        hms = hms-1;
    end

    % calculate mid stance of left front
    if X(16)<X(17)
        lf_ms = (X(16)+X(17))/2/X(22);
    elseif X(16)>X(17)
        lf_ms = (X(16)+X(17)+X(22))/2/X(22);
    end

    % calculate mid stance of right front
    if X(20)<X(21)
        rf_ms = (X(20)+X(21))/2/X(22);
    elseif X(20)>X(21)
        rf_ms = (X(20)+X(21)+X(22))/2/X(22);
    end

    % calculate midstance of hind leg pair
    if norm(lf_ms-rf_ms)<0.5
        fms = (lf_ms+rf_ms)/2;
    else
        fms = (lf_ms+rf_ms+1)/2;
    end

    if fms>1
        fms = fms-1;
    end

    MidstanceDiff = norm(fms-hms);


end