%%% This script measures the rupture length leading up to a step-over
clear all; close all;
%% load propagation length shapefiles 
stepover_prop_breached = shaperead('breached_stepover_leadup.shp');
stepover_prop_unbreached = shaperead('unbreached_stepover_leadup.shp');

%% load FDHI data with info 
FDHI_data = readtable('data_FDHI.xlsx');
event_info = readtable('event_info.csv');
event = FDHI_data.eq_name;
EQname_all = event_info.EQ_name;

%% load all step-over shapefiles

% Get a list of files in the current directory
files = dir('stepover*.shp');
Lmin = [];
width = [];
type = {}; 
style = {}; 
eventp = {};
label = {};


% Iterate over the shapefiles
for i = 1:length(files)
    % Get the filename
    filename = files(i).name;
    components = strsplit(filename,{'_','.'});
    % breached vs unbreached 
    shapefile_subtype = components{2};
    
   % breached
    % Load the shapefile
    shapefile = fullfile(pwd, filename);
    shapes = shaperead(shapefile);

    if isempty(shapes)
        continue
    % Process the shapefile as desired
    else
    for n=1:numel(shapes)
        xfault = shapes(n).X;
        yfault = shapes(n).Y;
        xfault = xfault(~isnan(xfault));
        yfault = yfault(~isnan(yfault));

    % get event name
    eventi = components{3};
    if length(components) > 4
    eventi = append(components{3},'_',components{4});
    else 
    end

    % get event coordinate system
    idx_event = find(strcmp(event,eventi));
    datasub = FDHI_data(idx_event,:);

    event_infosubidx = find(strcmp(EQname_all,eventi));
    event_infosub = event_info(event_infosubidx,:);
    
    eventp = [eventp; eventi];
    zone = datasub.zone;
    zone = zone{1};

if length(zone) == 3
    zone = cellstr(zone')';
    zone_n = append(zone{1},zone{2});
    zone_n = str2double(zone_n); 
    hem = zone{3};
    
elseif length(zone) == 2
    zone = cellstr(zone')';
    zone_n = str2double(zone{1}); 
    hem = zone{2};
else
    error('Length of zone string must be 2 or 3 characters')
end

[xfault,yfault]=wgs2utm(yfault,xfault,zone_n,hem);

stylei = event_infosub.style; 
style = [style; stylei];

labeli = shapefile_subtype; 
label = [label; labeli];

if isfield(shapes,'type')
    mech_type = shapes.type; 
    if strcmp(mech_type,'''T''')
     mech_typei = 'releasing';
    elseif strcmp(mech_type,'''C''')
     mech_typei = 'restraining';
    else 
    end
else
end

type = [type; mech_typei];

    if strcmp(shapefile_subtype,'breached')
         dist_pts =[];
         ikeep =[];
         iidx = [];

         Xb = []; 
         Yb = [];

    % find the two closest locations in the propagation lengths
    

    for p=1:numel(stepover_prop_breached)

        % coordinate transformation
        X_breached = stepover_prop_breached(p).X;
        X_breached = X_breached(~isnan(X_breached));
        Y_breached = stepover_prop_breached(p).Y;
        Y_breached = Y_breached(~isnan(Y_breached));
        [X_breached,Y_breached]=wgs2utm(Y_breached,X_breached,zone_n,hem);

        % find nearest neighbor in segment x and y to the step-over x and y
        [idx_pts,dist_ptsi] = dsearchn([xfault' yfault'],[X_breached' Y_breached']);
        dist_pts = [dist_pts; dist_ptsi];
        Xb = [Xb; X_breached'];
        Yb = [Yb; Y_breached'];
        ikeep = [ikeep; repelem(p,length(dist_ptsi))'];
    end 

    [sorted_vector, sorted_indices] = sort(dist_pts);
    min_first = sorted_indices(1);
    min_second = sorted_indices(2);

    dist_mina = sorted_vector(1);
    dist_minb = sorted_vector(2);

    if dist_mina>8000
        disp(eventi)
        error('Distance between point and step-over exceeds 10 km - a')
    elseif dist_minb>8000
        disp(eventi)
        error('Distance between point and step-over exceeds 10 km - b')
    else
    end

    Xkeep1 = Xb(min_first); % contains x coordinates of two nearest points. i.e. vertices of the two lines to measure
    Ykeep1 = Yb(min_first);   
    Xkeep2 = Xb(min_second); % contains x coordinates of two nearest points. i.e. vertices of the two lines to measure
    Ykeep2 = Yb(min_second);   
    ikeep1 = ikeep(min_first);
    ikeep2 = ikeep(min_second);


    if ikeep1==ikeep2
        disp(eventi)
        error('Two closest indices belong to same segment')
    
    else
       
    end

    % measure length of the two segments
    L1 = measure_length(stepover_prop_breached(ikeep1).X,stepover_prop_breached(ikeep1).Y,zone_n,hem);
    L2 = measure_length(stepover_prop_breached(ikeep2).X,stepover_prop_breached(ikeep2).Y,zone_n,hem);
    combined_L = [L1;L2];
    minL = min(combined_L);

    if L1 == L2
        error('Lenght of both segments the same, are you measuring the right thing?')
    else
    end 
    
    Lmin = [Lmin; minL];

    widthi = measure_length(xfault,yfault,zone_n,hem);
    width = [width; widthi];

    elseif strcmp(shapefile_subtype,'unbreached')
         dist_pts =[];
         ikeep =[];
         iidx = [];

         Xb = []; 
         Yb = [];

         
    % find the closest location in the propagation lengths
 for p=1:numel(stepover_prop_unbreached)
        % coordinate transformation
        X_unbreached = stepover_prop_unbreached(p).X;
        X_unbreached = X_unbreached(~isnan(X_unbreached));
        Y_unbreached = stepover_prop_unbreached(p).Y;
        Y_unbreached = Y_unbreached(~isnan(Y_unbreached));
        [X_unbreached,Y_unbreached]=wgs2utm(Y_unbreached,X_unbreached,zone_n,hem);

        % find nearest neighbor in segment x and y to the step-over x and y
        [idx_pts,dist_ptsi] = dsearchn([xfault' yfault'],[X_unbreached' Y_unbreached']);
        dist_pts = [dist_pts; dist_ptsi];
        Xb = [Xb; X_unbreached'];
        Yb = [Yb; Y_unbreached'];
        ikeep = [ikeep; repelem(p,length(dist_ptsi))'];
    end 

    [sorted_vector, sorted_indices] = sort(dist_pts);
    min_first = sorted_indices(1);
   
    Xkeep1 = Xb(min_first); % contains x coordinates of two nearest points. i.e. vertices of the two lines to measure
    Ykeep1 = Yb(min_first);   
    ikeep1 = ikeep(min_first);

    dist_min = sorted_vector(1);
    if dist_min>10000
        disp(eventi)
        scatter(Xkeep1,Ykeep1,'filled')
        hold on
        plot(xfault,yfault)
        error('Distance between point and step-over exceeds 10 km')
    else
    end

    % measure length of the two segments
    L1 = measure_length(stepover_prop_unbreached(ikeep1).X,stepover_prop_unbreached(ikeep1).Y,zone_n,hem);
    minL = L1;

    Lmin = [Lmin; minL];

    widthi = measure_length(xfault,yfault,zone_n,hem);
    width = [width; widthi];
    

    else
    error('Features must be breached or unbreached')
    end

    end

end
end

% export data
T = table(width, Lmin, type, style, label,eventp);
T.Properties.VariableNames = {'width', 'length', 'type', 'style', 'label','event'};
writetable(T, 'stepover_info.csv');

%% function dumpster

function [L] = measure_length(fault_x,fault_y,zone,hem)
fault_x = fault_x(~isnan(fault_x)); % removes NaN artifact at end of each fault in shapefile
fault_y =fault_y(~isnan(fault_y));
if fault_y<90
[fault_x,fault_y]=wgs2utm(fault_y,fault_x,zone,hem);
else
end
 % calculate length
x_1 = fault_x(1:end-1);
x_2 = fault_x(2:end);
y_1 = fault_y(1:end-1);
y_2 = fault_y(2:end);
segment_length = sqrt((x_1-x_2).^2+(y_1-y_2).^2); % note transformation to local coordinate system 
L = sum(segment_length);
end