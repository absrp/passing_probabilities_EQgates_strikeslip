%% this script takes in shapefiles of earthquake gates along surface ruptures in the FDHI
%% database and measures their geometrical attributes

% required inputs

% earthquake gate map shapefiles per event
% shapefile of ECS lines from FDHI database '_FDHI_FLATFILE_ECS_rev2.shp'
% info from FDHI appendix 'data_FDHI.xlsx'

% required functions (available from Mathworks in links below): 
% function wsg2utm (version 2) https://www.mathzworks.com/matlabcentral/fileexchange/14804-wgs2utm-version-2
% function distance2curve https://www.mathworks.com/matlabcentral/fileexchange/34869-distance2curve#:~:text=Distance2curve%20allows%20you%20to%20specify,and%20the%20closest%20point%20identified.

clear; close all;

% import FDHI data
FDHI_data = readtable('data_FDHI.xlsx');

% import shapefiles and extract information from all of them 
shapefiles = dir('*.shp'); % access all shapefile names in the folder

% create results table
all_results = table();

fault_x = [];
fault_y = [];

reflines_all = shaperead('_FDHI_FLATFILE_ECS_rev2.shp'); 

%% populate spreadsheet with geometry and other information

for i=30%:length(shapefiles)
% read shapefile
shapename = shapefiles(i).name;
maplines = shaperead(shapename); 

if isempty(maplines) 
    continue % skip for loop iteration if shapefile is empty
else 

name = strsplit(shapename,{'_','.'}); % string containing shapefile name

% extract info from shapefile name 

% feature type
shapefile_type = name{1};

% breached vs unbreached 
shapefile_subtype = name{2};
if strcmp(shapefile_subtype,'breached')
elseif strcmp(shapefile_subtype,'unbreached')
else
    error('Features must be breached or unbreached')
end

% earthquake name 
EQ_name= name{3};

if length(name) > 4
    EQ_name = append(name{3},'_',name{4});
else 
end

% find data associated with select earthquake
EQ_select = find(strcmp(FDHI_data.eq_name,EQ_name));
EQ_ID = FDHI_data.EQ_ID(EQ_select);
data = FDHI_data(EQ_select,:);
epicenter_xall = data.hypocenter_longitude_degrees;
epicenter_yall = data.hypocenter_latitude_degrees;
coordsx = data.longitude_degrees;
coordsy = data.latitude_degrees;
slip = data.recommended_net_preferred_for_analysis_meters;
epicenter_x = epicenter_xall(1);
epicenter_y = epicenter_yall(1);

% find ECS line for select earthquake
celllines = struct2cell(reflines_all)'; 
reflinesloc = find(cell2mat(celllines(:,5)) == EQ_ID(1)); 
reflines = reflines_all(reflinesloc);

%% measure length of maplines (i.e. gap and step-over length) from shapefile)
L_line = []; % create vector to store geometry data
measurement_type_line = {}; 
data = FDHI_data(EQ_select,:);

% extract utm zone from FDHI database
zone = data.zone;
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

% measure step-over width, gap length, and bend and splay angle
for n = 1:length(maplines) 
  [L_line(n),measurement_type_line{n}] =  measure_length_angle(maplines(n).X,maplines(n).Y,zone_n,hem,shapefile_type); 
  [fault_xi,fault_yi]= savecoords(maplines(n).X,maplines(n).Y,zone_n,hem);
  fault_x = [fault_x; fault_xi'];
  fault_y = [fault_y; fault_yi'];
end 

% save files as transtensional or compressional (based on identifier in
% shapefile)

mech_type = {};

for n = 1:length(maplines) 
    if isfield(maplines(n),'type')
        mech_type{n} = maplines(n).type; 
    else 
        mech_type{n} = 'NaN';
    end 
end

for p = 1:length(mech_type)
    if strcmp(mech_type{p},'''T''')
     mech_type{p} = 'releasing';
    elseif strcmp(mech_type{p},'''C''')
     mech_type{p} = 'restraining';
    else 
    end
end 

distance = [];
type_bend = {};
spacing = [];

% save double bend length and step-over proxy width

for n = 1:length(maplines) 
    if isfield(maplines(n),'distance')
        distance(n) = maplines(n).distance; 
        type_bend{n} = 'NaN'; 
        spacing(n) = 0;
        
   elseif strcmp(shapefile_type,'bend')
   for n = 1:length(maplines) 
       
    [distance(n),type_bend{n}] =  measure_length(maplines(n).X,maplines(n).Y,zone_n,hem);
    spacing(n) = measure_length_stepover_bend(maplines(n).X,maplines(n).Y,zone_n,hem);
    
    if strcmp(type_bend{n},'single')
        mech_type{n} = 'NaN';
        spacing(n) = 0;
    else 
    end 
    
   end       
    
    else
        distance(n) = 0;
        type_bend{n} = 'NaN'; 
        spacing(n) = 0;
    end 
end

xcheck = [];
ycheck = [];
latcheck = [];
loncheck = [];

for n=1:length(maplines)
    [xchecki, ychecki] = wgs2utm(maplines(n).Y,maplines(n).X,zone_n,hem);
    xcheck(n) = xchecki(1); 
    ycheck(n) = ychecki(1);
    latcheck(n) = maplines(n).Y(1); 
    loncheck(n) = maplines(n).X(1); 
end 

dimcheck = size(maplines);
dimcheck = dimcheck(:,1); 

%% location of gate along rupture (referenced to ECS files in FDHI database)

loc_along = [];
normalized_loc_along = [];
total_rupturelength = [];
distance_to_slipmax = [];

for n = 1:length(maplines)
   [total_rupturelengthi,loc_alongi,normalized_loc_alongi] = measure_location_along_rupture(maplines(n).X,maplines(n).Y,reflines.X,reflines.Y,zone_n,hem);
    loc_along = [loc_along; loc_alongi];
    normalized_loc_along = [normalized_loc_along; normalized_loc_alongi];
    total_rupturelength = [total_rupturelength; total_rupturelengthi];
end 

distance_to_epicenter = [];

for n = 1:length(maplines)
   [distance_to_epicenteri] = measure_distance_to_epicenter(maplines(n).X,maplines(n).Y,epicenter_x,epicenter_y,zone_n,hem);
    distance_to_epicenter = [distance_to_epicenter; distance_to_epicenteri];
end 

for n = 1:length(maplines)
    [distance_to_slipmaxi] = measure_distance_to_slipmax(maplines(n).X,maplines(n).Y,coordsx,coordsy,slip,zone_n,hem);
    distance_to_slipmax = [distance_to_slipmax; distance_to_slipmaxi];
end 

slip_at_gate = [];
normalized_slip_at_gate = [];

figure
for n = 1:length(maplines)
    [slip_at_gatei,normalized_slip_at_gatei] = find_slip_at_gate(maplines(n).X,maplines(n).Y,coordsx,coordsy,slip,zone_n,hem,shapefile_type);
    slip_at_gate = [slip_at_gate; slip_at_gatei];
    normalized_slip_at_gate = [normalized_slip_at_gate; normalized_slip_at_gatei];
end


%% extract info from the nearest data point near the earthquake gate from the FDHI database
% subset section of the FDHI database for desired earthquake

SRL_data = readtable('cumulative_displacements.xlsx'); 
eventSRL = SRL_data.Event;

idxSRL = find(strcmp(eventSRL,EQ_name));
if isempty(idxSRL)
    error('Earthquake name not in SRL database')
else 
end

Cumdisp = SRL_data.CumulativeDisplacement_km_(idxSRL); 

% slip = data.fps_central_meters;
magnitude = data.magnitude;
% fault_zone_width = data.fzw_central_meters;
% lithology = data.geology;
% coordsx = data.longitude_degrees;
% coordsy = data.latitude_degrees; 
date = data.eq_date;
hypo_lat = data.hypocenter_latitude_degrees;
hypo_lon = data.hypocenter_longitude_degrees;
EQ_style = data.style;
zone = data.zone;
zone = zone{1};

%% write data to table

allresults_i = table(...
    repelem(EQ_ID(1),dimcheck)', ...
    repelem(string(EQ_name),dimcheck)',...
    repelem(date(1),dimcheck)', ...
    repelem(magnitude(1),dimcheck)', ...
    repelem(EQ_style(1),dimcheck)',...
    repelem(hypo_lat(1),dimcheck)',...
    repelem(hypo_lon(1),dimcheck)',...
    repelem(string(shapefile_type),dimcheck)',...
    repelem(string(shapefile_subtype),dimcheck)',...
    string(mech_type)',...
    string(type_bend)',...
    distance',...
    L_line',...
    spacing',...
    measurement_type_line',...
    loc_along,...
    total_rupturelength,...
    normalized_loc_along,...
    distance_to_epicenter,...
    slip_at_gate,...
    normalized_slip_at_gate,...
    repelem(string(zone),dimcheck)', ...
    repelem(Cumdisp(1),dimcheck)');

all_results = [all_results; allresults_i];

disp(EQ_name); % keeps track of progress
title(EQ_name)
end
end


%% export results

% assign header to table
all_results.Properties.VariableNames = {'FDHI ID',...
    'Earthquake',...
    'Date',...
    'Magnitude',...
    'Style',...
    'Hypocenter lat',...
    'Hypocenter lon',...
    'Feature',...
    'Breached or unbreached',...
    'Type (releasing or restraining)',...
    'Type (single or double)',...
    'Distance splay or double bend (m)',...
    'Length (m) or angle (deg)',...
    'Spacing double bend (m)',...
    'Type (length or angle)',...
    'Location along rupture',...
    'Total rupture length',...
    'Normalized location',...
    'Distance to epicenter',...
    'Slip at gate (m)',...
    'Normalized slip at gate',...
    'UTM zone',...
    'Cumulative displacement'};

% export file as csv 
writetable(all_results,'aEQgate_geometries.csv'); 

%% function dumpster
% functions that are called in the script go here 
function [fault_x,fault_y] = savecoords(fault_x,fault_y,zone,hem)
fault_x = fault_x(~isnan(fault_x)); % removes NaN artifact at end of each fault in shapefile
fault_y =fault_y(~isnan(fault_y));
if fault_y<90
[fault_x,fault_y]=wgs2utm(fault_y,fault_x,zone,hem);
else
end
end
function [L,measurement_type_line] = measure_length_angle(fault_x,fault_y,zone,hem,shapefile_type)
fault_x = fault_x(~isnan(fault_x)); % removes NaN artifact at end of each fault in shapefile
fault_y =fault_y(~isnan(fault_y));
if fault_y<90
[fault_x,fault_y]=wgs2utm(fault_y,fault_x,zone,hem);
else
end

% measure angle or length depending on shapefile type 

if strcmp(shapefile_type,'splay') % check if shapefile type is a splay
% measure angle between maplines
v1=[fault_x(1),fault_y(1)]-[fault_x(2),fault_y(2)];
v2=[fault_x(end),fault_y(end)]-[fault_x(2),fault_y(2)];
L=acos(sum(v1.*v2)/(norm(v1)*norm(v2)));
L = rad2deg(L);
measurement_type_line = 'angle';

elseif strcmp(shapefile_type,'bend') % check if shapefile type is a bend
    if length(fault_x) == 3
    % measure angle between maplines
        v1=[fault_x(2),fault_y(2)]-[fault_x(1),fault_y(1)];
        v2=[fault_x(end),fault_y(end)]-[fault_x(2),fault_y(2)];
        L=acos(sum(v1.*v2)/(norm(v1)*norm(v2)));
        L = rad2deg(L);
        measurement_type_line = 'angle'; % single bend
        if L>90
            L = 180-L;
        elseif L>180
            error('Angle larger than 180')
        elseif L<0
            error('Angle smaller than 0')
        else
            L = L;
        end
    elseif length(fault_x) == 4
        v1=[fault_x(2),fault_y(2)]-[fault_x(1),fault_y(1)];
        v2=[fault_x(3),fault_y(3)]-[fault_x(2),fault_y(2)];
        anga=acos(sum(v1.*v2)/(norm(v1)*norm(v2)));
        anga = rad2deg(anga);
        measurement_type_line = 'angle'; % double bend

        % angle b test for double bends - ensuring the two angles in the double bend are not
        % too far apart from each other
        v1=[fault_x(3),fault_y(3)]-[fault_x(2),fault_y(2)];
        v2=[fault_x(4),fault_y(4)]-[fault_x(3),fault_y(3)];
        angb=acos(sum(v1.*v2)/(norm(v1)*norm(v2)));
        angb = rad2deg(angb);
        if anga-angb>10
            disp(anga-angb)
            L = (anga+angb)/2;
        else
            L = (anga+angb)/2;
        end
        
    else 
        disp(length(fault_x))
    error('Bends must contain three or four x,y coordinate pairs')

    end
        
elseif strcmp(shapefile_type,'stepover') % check if shapefile type is a step-over
% calculate length
x_1 = fault_x(1:end-1);
x_2 = fault_x(2:end);
y_1 = fault_y(1:end-1);
y_2 = fault_y(2:end);
segment_length = sqrt((x_1-x_2).^2+(y_1-y_2).^2); % note transformation to local coordinate system 
L = sum(segment_length);
measurement_type_line = 'length';

elseif strcmp(shapefile_type,'strand') % check if shapefile type is a step-over
% calculate length
x_1 = fault_x(1:end-1);
x_2 = fault_x(2:end);
y_1 = fault_y(1:end-1);
y_2 = fault_y(2:end);
segment_length = sqrt((x_1-x_2).^2+(y_1-y_2).^2); % note transformation to local coordinate system 
L = sum(segment_length);
measurement_type_line = 'length';

elseif strcmp(shapefile_type,'gap')  % check if shapefile type is a gap
% calculate length
x_1 = fault_x(1:end-1);
x_2 = fault_x(2:end);
y_1 = fault_y(1:end-1);
y_2 = fault_y(2:end);
segment_length = sqrt((x_1-x_2).^2+(y_1-y_2).^2); % note transformation to local coordinate system 
L = sum(segment_length);
measurement_type_line = 'length';

else
    error('ERROR: Shapefile type must be splay, gap, bend, or step-over (stepover)')
end 
end 
function [distance,bend_type] = measure_length(fault_x,fault_y,zone,hem)
fault_x = fault_x(~isnan(fault_x)); % removes NaN artifact at end of each fault in shapefile
fault_y =fault_y(~isnan(fault_y));
if fault_y<90
[fault_x,fault_y]=wgs2utm(fault_y,fault_x,zone,hem);
else
end
if length(fault_x) == 4
        segment_length = sqrt((fault_x(2)-fault_x(3)).^2+(fault_y(2)-fault_y(3)).^2); % note transformation to local coordinate system 
        distance = sum(segment_length);
        bend_type = 'double';
else
    distance = 0;
    bend_type = 'single';
end

end 
function [spacing] = measure_length_stepover_bend(fault_x,fault_y,zone,hem)
fault_x = fault_x(~isnan(fault_x)); % removes NaN artifact at end of each fault in shapefile
fault_y =fault_y(~isnan(fault_y));

if fault_y<90
[fault_x,fault_y]=wgs2utm(fault_y,fault_x,zone,hem);
else
end

if length(fault_x) == 4
    
        v1=[fault_x(2),fault_y(2)]-[fault_x(1),fault_y(1)];
        v2=[fault_x(3),fault_y(3)]-[fault_x(2),fault_y(2)];
        angle_rad = acos(sum(v1.*v2)/(norm(v1)*norm(v2)));
        angle = rad2deg(angle_rad);
        
        segment_length = sqrt((fault_x(2)-fault_x(3)).^2+(fault_y(2)-fault_y(3)).^2); % note transformation to local coordinate system 
        hypothenuse = sum(segment_length);
        
        % calculate step-over spacing for double bend
        spacing = sind(angle)*hypothenuse; 
        
        
else
    spacing = 0;
end

end 
function [total_rupturelength,loc_along,normalized_loc_along] = measure_location_along_rupture(fault_x,fault_y,refline_x,refline_y,zone,hem)

fault_x = fault_x(~isnan(fault_x)); % removes NaN artifact at end of each fault in shapefile
fault_y = fault_y(~isnan(fault_y));
refline_x = refline_x(~isnan(refline_x));
refline_y = refline_y(~isnan(refline_y));

[coords_gatex, coordsgatey] = wgs2utm(fault_y(1),fault_x(1),zone,hem); % first point on the gate
coords_gate = [coords_gatex' coordsgatey'];

[curvexy_x, curvexy_y] = wgs2utm(refline_y,refline_x,zone,hem);

% total length
x_1 = curvexy_x(1:end-1);
x_2 = curvexy_x(2:end);
y_1 = curvexy_y(1:end-1);
y_2 = curvexy_y(2:end);
segment = sqrt((x_1-x_2).^2+(y_1-y_2).^2); % note transformation to local coordinate system 
total_rupturelength = sum(segment);

spacing = 10; % discretizing rupture into 10 m spaced increments to resample
pt = interparc(0:(spacing/total_rupturelength):1,curvexy_x,curvexy_y,'linear'); 
pt_x = pt(:,1);
pt_y = pt(:,2);
curvexy_dense = [pt_x pt_y];

% a few tests to make sure everything is running ok
% disp('Number of points in ECS line')
% disp(length(refline_x))
% disp('Number of points in interpolated line')
% disp(length(pt_y))
% disp('Surface rupture length')
% disp(total_rupturelength)

[xy,dist,~] = distance2curve(curvexy_dense,coords_gate,'spline'); % find minimum distance between gate and ECS trace
[locpt,~] = dsearchn(curvexy_dense,xy); % finding what index corresponds to the location of the gate

% a few tests to make sure everything is running ok
% disp('Location on ECS trace that is closest to gate')
% disp(xy)
% disp('Chosen coordinates based on densified ECS curve - x')
% disp(pt_x(locpt))
% disp('Chosen coordinates based on densified ECS curve - y')
% disp(pt_y(locpt))
% disp('Distance between the gate and chosen coordinate')
% disp(dist)

% segment length
x_1 = pt_x(1:locpt-1);
x_2 = pt_x(2:locpt);
y_1 = pt_y(1:locpt-1);
y_2 =  pt_y(2:locpt);
segment = sqrt((x_1-x_2).^2+(y_1-y_2).^2); 
loc_along= sum(segment);

normalized_loc_along = loc_along/total_rupturelength; 
end
function [distance_to_epi] = measure_distance_to_epicenter(fault_x,fault_y,epi_x,epi_y,zone,hem)

fault_x = fault_x(~isnan(fault_x)); % removes NaN artifact at end of each fault in shapefile
fault_y = fault_y(~isnan(fault_y));

% hold on
% scatter(fault_x,fault_y,'b')
% scatter(epi_x,epi_y,'r')

[coords_gatex, coordsgatey] = wgs2utm(fault_y(1),fault_x(1),zone,hem);
coords_gate = [coords_gatex' coordsgatey'];

[hypoxy_x, hypoxy_y] = wgs2utm(epi_y,epi_x,zone,hem);
hypoxy = [hypoxy_x' hypoxy_y'];

[~,distance_to_epi] = dsearchn(hypoxy,coords_gate); % find minimum distance between gate and epicenter

end 
function [distance_to_slipmax] = measure_distance_to_slipmax(fault_x,fault_y,coords_x,coords_y,slip,zone,hem)

% find location of maximum displacement
[coordsx_slip,coordsy_slip] = wgs2utm(coords_y,coords_x,zone,hem);

fault_x = fault_x(~isnan(fault_x)); % removes NaN artifact at end of each fault in shapefile
fault_y = fault_y(~isnan(fault_y));
[coords_gatex, coordsgatey] = wgs2utm(fault_y(1),fault_x(1),zone,hem);
idxmax = find(slip==max(slip));
if idxmax>2
    idxmax = idxmax(1);
end 

locx_maxslip = coordsx_slip(idxmax);
locy_maxslip = coordsy_slip(idxmax);

% measure distance between gate and coordinates of maximum slip
coordsgate = [coords_gatex coordsgatey];
coordsslip = [locx_maxslip locy_maxslip];
%coords_measure = [coordsgate; coordsslip];
[~,distance_to_slipmax] = dsearchn(coordsgate,coordsslip);
%distance_to_slipmax = pdist(coords_measure);

end 
function [slip_at_gate,normalized_slip_at_gate] = find_slip_at_gate(fault_x,fault_y,coords_x,coords_y,slip,zone,hem,feature)
% find location of all slip measurements
[coordsx_slip,coordsy_slip] = wgs2utm(coords_y,coords_x,zone,hem);
slip_nonzero =  find(slip>-0.1); % avoid -999 errors while including measures of zero slip
coordsx_slip = coordsx_slip(slip_nonzero);
coordsy_slip = coordsy_slip(slip_nonzero);
slip = slip(slip_nonzero);

% find location of gate
fault_x = fault_x(~isnan(fault_x)); % removes NaN artifact at end of each fault in shapefile
fault_y = fault_y(~isnan(fault_y));
[coords_gatex, coords_gatey] = wgs2utm(fault_y,fault_x,zone,hem);

% measure distance between gate and each slip point 
coordsgate = [coords_gatex' coords_gatey'];
coordsslip = [coordsx_slip coordsy_slip];

radius = 500; % meters

% from middle point of bend 
if strcmp(feature,'bend')
    if length(coords_gatex) == 3 % single bend
    [distance_to_slip] = pdist2(coordsgate(2,:),coordsslip);
    idx_radius = find(distance_to_slip<radius);
    length(idx_radius)
    slip_in_radius = slip(idx_radius);


    elseif length(coords_gatex) == 4 % double bend
    %[distance_to_slip] = pdist2(coordsgate(2:3,:),coordsslip);
    [distance_to_slip_pt2, I_pt2] = pdist2(coordsgate,coordsslip(2,:),'euclidean', 'smallest', 2); % for bend length pt 1
    [distance_to_slip_pt3, I_pt3] = pdist2(coordsgate,coordsslip(3,:),'euclidean', 'smallest', 2); % for bend length pt 2

    idx_radius_pt2 = find(distance_to_slip_pt2<radius);
    idx_radius_pt3 = find(distance_to_slip_pt3<radius);
    I_pt2 = I_pt2(idx_radius_pt2);
    I_pt3 = I_pt3(idx_radius_pt3);
    % find overlap between I_pt2 and I_pt3
    commonValues = intersect(I_pt2, I_pt3);
    [~, loc2] = ismember(I_pt2, commonValues);
    I_pt2(loc2 ~= 0) = [];  
    [~, loc3] = ismember(I_pt3, commonValues);
    I_pt3(loc3 ~= 0) = []; 
    slipvals_idx = [I_pt2;I_pt3];

    slip_in_radius = slip(slipvals_idx);
    hold on
    scatter(coordsgate(2:3,1),coordsgate(2:3,2),'filled','MarkerFaceAlpha',0.1)
    scatter(coordsslip(:,1),coordsslip(:,2),'filled','k')
    scatter(coordsslip(slipvals_idx,1),coordsslip(slipvals_idx,2),'filled','m')


    else 
        error('Length of bend vector must be 3 or 4 elements')
    end

else % all other features

%[distance_to_slip] = pdist2(coordsgate,coordsslip);
    [distance_to_slip_pt2, I_pt2] = pdist2(coordsgate,coordsslip(2,:),'euclidean', 'smallest', 2); 
    [distance_to_slip_pt1, I_pt1] = pdist2(coordsgate,coordsslip(1,:),'euclidean', 'smallest', 2); 
    idx_radius_pt2 = find(distance_to_slip_pt2<radius);
    idx_radius_pt1 = find(distance_to_slip_pt1<radius);
    I_pt2 = I_pt2(idx_radius_pt2);
    I_pt1 = I_pt1(idx_radius_pt1);
    % find overlap between I_pt2 and I_pt3
    commonValues = intersect(I_pt2, I_pt1);
    [~, loc2] = ismember(I_pt2, commonValues);
    I_pt2(loc2 ~= 0) = [];  
    [~, loc1] = ismember(I_pt1, commonValues);
    I_pt1(loc1 ~= 0) = []; 
    slipvals_idx = [I_pt2;I_pt1]

    slip_in_radius = slip(slipvals_idx);
    hold on
    scatter(coordsgate(:,1),coordsgate(:,2),'filled','MarkerFaceAlpha',0.1)
    scatter(coordsslip(:,1),coordsslip(:,2),'filled','k')
    scatter(coordsslip(slipvals_idx,1),coordsslip(slipvals_idx,2),'filled','m')


end 

if length(slip_in_radius)>1
    slip_at_gate = mean(slip_in_radius); 
    % normalize slip at gate
    max_slip = max(slip); % max slip for event gate belongs to
    normalized_slip_at_gate = slip_at_gate/max_slip; 

else
    slip_at_gate = NaN;
    normalized_slip_at_gate = NaN;
end 



end