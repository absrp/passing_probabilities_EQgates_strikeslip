%% this script calculates the distance between the centroid of each zone of geometrical complexity and the distance to its nearest neighbor, examining the distribution of distances
% required inputs

% geometrical complexity shapefiles per event 
% info from FDHI appendix 'data_FDHI.xlsx' (from Sarmiento et al., 2021)

%% set-up 
clear all; close all 

currentDir = pwd;
addpath(genpath(fullfile(currentDir, 'Source_code')));
addpath(genpath(fullfile(currentDir, '11095762'))); % Zenodo repo data -- see required inputs above 

shapefileDir = fullfile(currentDir, '11095762/primary_EQgate_shapefiles_v1'); % shapefile directory in folder 11095762

%% load data

% import shapefiles
shapefiles = dir(fullfile(shapefileDir, '*.shp'));

shapefile_names = {};
shapefile_type = {};
shapefile_BU = {};
shapefile_event = {};

for n=1:numel(shapefiles)
    shapefile_namesi = shapefiles(n).name; 
    shapefile_names{n} = shapefile_namesi;
    
    % break down into string components
    name = strsplit(shapefile_namesi,{'_','.'}); % string containing shapefile name

    shapefile_type{n} = name{1};
    shapefile_BU{n} = name{2};
    shapefile_event{n} = name{3};
    
end

% event info from FDHI database
FDHI_data = readtable('data_FDHI.xlsx');
EQ = FDHI_data.eq_name;
[EQ,iEQ] = unique(EQ,'legacy'); 
type = FDHI_data.style(iEQ,:);
zone = FDHI_data.zone(iEQ,:);

locs = find(strcmp(type,'Strike-Slip'));
names = EQ(locs);
zone = zone(locs);


%% find centroids for each earthquake gate 

% for populating centroid table 
centroidx = [];
centroidy = [];
event_name = [];
length_gate = []; % for trouble shooting

for i = 1:numel(names) 
    namei = names(i);
     
    % deal with utm zone
    zonei = zone(i);
    zonei = zonei{1};

if length(zonei) == 3
    zonei = cellstr(zonei')';
    zone_n = append(zonei{1},zonei{2});
    zone_n = str2double(zone_n); 
    hem = zonei{3};
    
elseif length(zonei) == 2
    zonei = cellstr(zonei')';
    zone_n = str2double(zonei{1}); 
    hem = zonei{2};
else
    error('Length of zone string must be 2 or 3 characters')
end
    
    % find shapefiles with same event name
    locname = find(strcmp(shapefile_event,namei));
    
    % create structure with shapefiles from same event
    lines= {};
    for n=1:length(locname)
        idn = locname(n);
        variablename = shapefiles(idn).name;
        nametest = strsplit(variablename,{'_','.'}); % string containing shapefile name
        if strcmp(nametest{2},'breached') && (strcmp(nametest{1}, 'bend') || strcmp(nametest{1}, 'stepover')) || strcmp(nametest{1}, 'gap')% choosing only EQ gates
        lines{n} = shaperead(variablename);      
        else
            continue
        end
    end
    
    % compile all the centroids in shapefile 

    for n=1:length(lines)
        selected = lines{n};
        for p=1:length(selected)
            fault_x = selected(p).X;
            fault_y = selected(p).Y;
            fault_x = fault_x(~isnan(fault_x)); % removes NaN artifact at end of each fault in shapefile
            fault_y = fault_y(~isnan(fault_y)); 
            [fault_x, fault_y] = wgs2utm(fault_y,fault_x,zone_n,hem);
            if length(fault_x) == 4 % double bend
                P1 = [fault_x(2), fault_y(2)];
                P2 = [fault_x(3), fault_y(3)];
                midpoint = (P1(:) + P2(:)).'/2;
                centroidxi = midpoint(1); 
                centroidyi = midpoint(2); 
            elseif length(fault_x) == 2 % step-over, strand, and gap
                P1 = [fault_x(1), fault_y(1)];
                P2 = [fault_x(2), fault_y(2)];
                midpoint = (P1(:) + P2(:)).'/2;
                centroidxi = midpoint(1); 
                centroidyi = midpoint(2);               
            else % splay and single bend
                centroidxi = fault_x(2); 
                centroidyi = fault_y(2); 
            end
            event_name_i = repmat(namei,length(centroidyi));


            centroidx = [centroidx; centroidxi'];
            centroidy = [centroidy; centroidyi'];    
            event_name = [event_name; event_name_i']; % repeat event name times the number of centroid for book keeping
            length_gate = [length_gate; length(fault_x)'];
        end

    end
end 


%% measure the distance between each earthquake gate and its nearest neighbor (centroids)

distances = []; % to save nearest neighbor distances
centroid_table = table(centroidx, centroidy, event_name, 'VariableNames', {'centroid_x', 'centroid_y', 'event_name'});

for i=1:length(EQ)
% find indices of events
selected_event = EQ(i);
event_rows = find(strcmp(event_name,selected_event));
centroidx_event = centroidx(event_rows);
centroidy_event = centroidy(event_rows);
length_rows = length_gate(event_rows);

for c=1:length(centroidx_event)
    idxrem = find(centroidx_event == centroidx_event(c) & centroidy_event == centroidy_event(c));
    centroidsx_subset = centroidx_event(1:end ~= idxrem);
    centroidsy_subset = centroidy_event(1:end ~= idxrem);
    [k,disti] = dsearchn([centroidsx_subset,centroidsy_subset],[centroidx_event(c), centroidy_event(c)]);
    if disti>10^4
        disp(selected_event) % catch problematic centroids and examine them
        disp(length_rows(c))
        disp(centroidx_event)
        disp(centroidy_event)
        disp(c)
    else
        
    end
    distances = [distances; disti'];
end 

end


%% visualize and fit distribution of nearest neighbors

    figure
    histogram(distances,20,'FaceColor',[0.9294    0.6941    0.1255])
    xlabel('Distance to nearest neighbor breached gate (m)')
    ylabel('Frequency')
    set(gca,'FontSize',14)

disp('Mean distance in meters:')
disp(mean(distances))

%% ECDF
figure()
[F, X] = ecdf(distances);
plot(X, 1-F,'Color','k','linewidth',2);
xlabel('Distance to nearest neighbor (m)');
ylabel('1 - Cumulative Probability');
set(gca,'FontSize',14,'XScale','log')

% CDFs (log normal, exponential, and Weibull)
pd_lognormal = fitdist(distances, 'LogNormal');
pd_exponential = fitdist(distances, 'Exponential');
pd_weibull = fitdist(distances, 'Weibull');
x = linspace(min(distances),max(distances),10000);
fittedCDF_lognormal = cdf(pd_lognormal, x);
fittedCDF_exponential = cdf(pd_exponential, x);
fittedCDF_weibull = cdf(pd_weibull, x);
hold on;
plot(x, 1-fittedCDF_lognormal, 'Color',[0.8510    0.3255    0.0980], 'LineWidth', 2);  
plot(x, 1-fittedCDF_exponential, 'Color',[0.4667    0.6745    0.1882], 'LineWidth', 2);  
plot(x, 1-fittedCDF_weibull, 'Color',[0.9294    0.6941    0.1255], 'LineWidth', 2);  

legend('Empirical CDF', 'log-normal', 'exponential', 'Weibull');
xlim([0,12000])

% CDF gamma
% params = gamfit(dist');
% gamma_cdf = gamcdf(x, params(1), params(2)); % Use the estimated parameters
% gamma_cdf_plot = plot(x, 1-gamma_cdf, 'Color', 'b', 'LineWidth', 2);

% saveas(gcf,'a_lognormal_exp_Weibull_CDF.pdf')