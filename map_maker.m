clear all; close all; 

secondary_fractures = shaperead('_FDHI_FLATFILE_RUPTURES_20220719.shp');
event_info = readtable('event_info.txt');
zone = event_info.zone; 
event = event_info.event;
region = event_info.region; 

rupturemaps = shaperead('_FDHI_FLATFILE_RUPTURES_20220719.shp');
Qfaults = shaperead('Qfaults_cut_fdhi.shp'); % cut Qfaults
GEM = shaperead('GEM_cut_FDHI.shp');   
AFEAD = shaperead('AFEAD_FDHI.shp');
NZ = shaperead('NZAFD_Oct_2020.shp');

% find event names in shapefile order
shapefiles = dir('*.shp'); % access all shapefile names in the folder

names = {};
for n=1:length(shapefiles)
namesi = shapefiles(n).name;
strpl = strsplit(namesi,{'_','.'});
names{n} = strpl{3};
end

 %%

for i= 1:length(event)
fig = figure
    eventi = event(i); 
    
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

% regional rupture map 
regioni = region(i); 

if strcmp(regioni,'Qfaults')
    lines_regional = Qfaults; 
elseif strcmp(regioni,'AFEAD')
    lines_regional = AFEAD;
elseif strcmp(regioni,'GEM')
    lines_regional = GEM;
else 
    lines_regional = NZ;
end 


for n=1:numel(lines_regional)
    plotline(lines_regional(n).X,lines_regional(n).Y,zone_n,hem,[0.5 0.5 0.5])
    hold on 
end

% find shapefiles with same event name
celllines = struct2cell(rupturemaps)'; 
reflinesloc = find(strcmp(celllines(:,6),eventi)); 
reflines = rupturemaps(reflinesloc);

% plot surface rupture map
for n=1:numel(reflines)
    plotline(reflines(n).X,reflines(n).Y,zone_n,hem,'k')
    hold on 
end

% x and y limits in agreement with event
xevent = [reflines(:).X];
yevent = [reflines(:).Y];
[xevent,yevent] = wgs2utm(yevent,xevent,zone_n,hem);


% plot earthquake gates over map 
  % find shapefiles with same event name
    idn = strcmp(names,event(i)); % find locations where name of event and name of shapefile overlap 
    idnx = find(idn == 1);
    
    for p=1:length(idnx)
        selectshp = idnx(p);
        variablename = shapefiles(selectshp).name;
        nametest = strsplit(variablename,{'_','.'}); % string containing shapefile name
        lines = shaperead(variablename);
        
    for n=1:length(lines)    
        if isempty(lines) 
            continue % skip for loop iteration if shapefile is empty
        else 
    
        if strcmp(nametest{1},'gap')
            if strcmp(nametest{2},'breached')
                plotline(lines(n).X,lines(n).Y,zone_n,hem,[0.9294    0.6941    0.1255])
            else 
                plotline(lines(n).X,lines(n).Y,zone_n,hem,[1     1     0])
            end               
        elseif strcmp(nametest{1},'bend')
            if strcmp(nametest{2},'breached')
                plotline(lines(n).X,lines(n).Y,zone_n,hem,[0.4667    0.6745    0.1882])
            else 
                plotline(lines(n).X,lines(n).Y,zone_n,hem,[0.0353    0.5098    0.1608])
            end    
        elseif strcmp(nametest{1},'stepover')
            if strcmp(nametest{2},'breached')
                plotline(lines(n).X,lines(n).Y,zone_n,hem,[0    0.4471    0.7412])
            else 
                plotline(lines(n).X,lines(n).Y,zone_n,hem,[0    0.8667    1.0000])
            end   
        elseif strcmp(nametest{1},'strand') 
            plotline(lines(n).X,lines(n).Y,zone_n,hem,[0.4941    0.1843    0.5569])
        elseif strcmp(nametest{1},'splay')
            if strcmp(nametest{2},'breached')
                plotline(lines(n).X,lines(n).Y,zone_n,hem,[1     0     0])
            else 
                plotline(lines(n).X,lines(n).Y,zone_n,hem,[0.6353    0.0784    0.1843])
            end   
        else
            continue
        end
        end
    end

set(gca,'YTickLabel',[],'XTickLabel',[])
title(event(i))
axis equal

 deltax = max(xevent)-min(xevent);
    deltay = max(yevent)-min(yevent);

    xlim([min(xevent)-deltax/2 max(xevent)+deltax/2]); 
    ylim([min(yevent)-deltay/2 max(yevent)+deltay/2]);  

scalebar

set(gcf,'PaperOrientation','landscape'); 
pdf_name = strcat(char(event(i)),'.pdf');
print (pdf_name, '-dpdf', '-fillpage', '-r300')

end
end
%% function dumpster

function [] = plotline(x,y,zone_n,hem,col)
[fault_x,fault_y]= wgs2utm(y,x,zone_n,hem);
plot(fault_x,fault_y,'Color',col)
axis equal
end 