% This script converts the kmz files in the FDHI database into shapefiles,
% keeping only the fields of interest. 
clear all; close all;
% required inputs
% kmz files per event, downloadable from the FDHI database appendices

% required functions (available from Mathworks in links below): 
% function kmz2struct
% (https://www.mathworks.com/matlabcentral/fileexchange/102609-kmz2struct/?s_tid=LandingPageTabfx)

% load all kmz files in current directory 
kmz_files = dir('*.kmz');
field = 'Line';
subfield = '/Ruptures - Principal';

fields_to_remove = {'Color', 'Description'};

%% loop through files 

for i=1:numel(kmz_files)

 % load kmz file and select fields 
 filename = kmz_files(i).name;
 struct_event = kmz2struct(filename);
 subset_indices = strcmp({struct_event.Geometry}, field);
 subset_struct = struct_event(subset_indices);
 event_name = subset_struct(1).Name;
 
 % account for typos in FDHI files and some events containing multiple
 % primary and secondary datasets, or different labels
 if strcmp(event_name, 'Duzce')
    subfield = '/Rupturs - Principal';
    subset_indices = strcmp({subset_struct.Folder}, subfield);
    subset_struct = subset_struct(subset_indices);
 elseif strcmp(event_name,'Darfield')
    subfield = '/Ruptures (RUP_DS_ID=103) - Principal';
    subset_indices = strcmp({subset_struct.Folder}, subfield);
    subset_struct = subset_struct(subset_indices);
 elseif strcmp(event_name,'Ridgecrest1') || strcmp(event_name,'Ridgecrest2')
    folders = {subset_struct.Folder};
    index_145 = strcmp(folders, '/Ruptures (RUP_DS_ID=145) - Principal');
    index_132 = strcmp(folders, '/Ruptures (RUP_DS_ID=132) - Principal');
    subset_indices = find(index_145 | index_132);
    subset_struct = subset_struct(subset_indices);
 elseif strcmp(event_name,'IzuPeninsula')
    subfield = '/Ruptues - Principal';
    subset_indices = strcmp({subset_struct.Folder}, subfield);
    subset_struct = subset_struct(subset_indices);    
 else 
    subfield =  '/Ruptures - Principal';
    subset_indices = strcmp({subset_struct.Folder}, subfield);
    subset_struct = subset_struct(subset_indices);
 end

 % remove fields not needed for shapefile 
 subset_struct = rmfield(subset_struct, fields_to_remove);

 % turn into shapefile 
 shp_extension = '_ruptures.shp';
 shapename = [event_name, shp_extension];
 shapewrite(subset_struct, shapename);

 % disp(i) % for debugging

end
