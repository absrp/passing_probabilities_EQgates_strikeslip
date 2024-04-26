function kmz_struct = kmz2struct(filename)
% Load a KML or a KMZ file as a MATLAB structure.
%
% :param filename: Filename of the KML or KMZ file to load.
% :param kmz_struct: Structure with data for the Polygon, Line and Point
% geometry placemarks found in the "filename."
    % Determine if filename is a kml or kmz file.
    [~,~,ext] = fileparts(filename);
    if strcmpi(ext,'.kmz')
        % Create temp .kmz2struct directory
        userDir = fullfile( ...
            char(java.lang.System.getProperty('user.home')), ...
            '.kmz2struct');
        if ~exist(userDir,'dir')
            mkdir(userDir);
        end
        % Upzip the file to the temp directory
        unzip(filename, userDir);
        % Search the temp directory for kml files and load each kml file
        % individually
        files = dir(fullfile(userDir, '**', '*.kml'));
        N = length(files);
        kmlStructs = cell([1 N]);
        for i = 1:length(files)
            kmlStructs{i} = readKMLfile( ...
                fullfile(files(i).folder, files(i).name));
        end
        kmz_struct = vertcat(kmlStructs{:});
        % Clean up the temp directory
        rmdir(userDir,'s');
    else
        % If input is a kml load that kml
        kmz_struct = readKMLfile(filename);
    end
end
function kmlStruct = readKMLfile(filename)
% Function to read in kml file.
%
% :param filename: KML filename to load.
    % Read the KML as an XML DOM
    doc = xmlread(filename);
    starting_node =  doc.item(0).item(1);
    % Find styles and create a hash table with a struct of colors for the
    % different geometries
    styles = starting_node.getElementsByTagName('Style');
    style_hash = containers.Map('KeyType','char','ValueType','Any');
    for j = 0:styles.getLength-1
        idxml = styles.item(j).getAttributes.getNamedItem('id');
        if ~isempty(idxml)
            id = char(idxml.getTextContent);
            [pointcolor,linecolor,polycolor] = parseStyle(styles.item(j));
            style_hash(id) = struct('pointcolor',pointcolor,...
                'linecolor',linecolor,...
                'polycolor',polycolor);
        end
    end
    % Handle style maps
    style_maps = starting_node.getElementsByTagName('StyleMap');
    for j = 0:style_maps.getLength-1
        id = char(style_maps.item(j).getAttributes.getNamedItem('id').getTextContent);
        keys = style_maps.item(j).getElementsByTagName('key');
        index = 0;
        for k = 0:keys.getLength-1
            found = strcmp(char(keys.item(k).getTextContent),'normal');
            if found; index = k; break; end
        end
        url_node =keys.item(index).getParentNode.getElementsByTagName('styleUrl');
        style_url = char(url_node.item(0).getTextContent);
        style_hash(id) = style_hash(style_url(2:end));
    end
    % Find Geometries
    kmlStruct = recursive_kml2struct(starting_node,'',style_hash);
end
function kmlStruct = recursive_kml2struct(folder_node,folder,style_hash)
% Recursive Function to search for kml placemarks
%
% :param folder_node: Node to search for placemarks in.
% :param folder: Character array to represent the folder of folder_node.
% :param style_hash: Hash table of structures containing possible colors
% for the kml placemarks.
    % Find number of placemarks and name of folder.
    name = 'none';
    number_placemarks = 0;
    for i = 0:folder_node.getLength()-1
        if strcmp(folder_node.item(i).getNodeName,'Placemark')
            number_placemarks = number_placemarks + 1;
        elseif strcmp(folder_node.item(i).getNodeName,'name')
            name = char(folder_node.item(i).getTextContent);
        end
    end
    % If we are searching a folder set the new folder name.
    if strcmpi(folder_node.getNodeName,'Folder')
        folder = [folder '/' name];
    end
    % Find Placemark Data
    count = 1;
    kmlStructs = cell([1 number_placemarks]);
    for i = 0:folder_node.getLength()-1
        current = folder_node.item(i);
        NodeName = current.getNodeName;
        if strcmpi(NodeName,'Folder')
            kmlStructs{count} = recursive_kml2struct(current,folder,style_hash);
            count = count + 1;
        elseif strcmpi(NodeName,'Placemark')
            kmlStructs{count} = parsePlacemark(current,folder,style_hash);
            count = count + 1;
        end
    end
    kmlStruct = horzcat(kmlStructs{:});
end
function kmlStruct = parsePlacemark(node,folder,style_hash)
% Parse a KML placemark DOM element to a MATLAB struct.
%
% :param node: The KML placemark node.
% :param folder: Character array reperesenting the folder of the placemark.
% :param style_hash: A hash table containing styles they may apply to this
% place mark.
    % Get the name of this element
    namexml = node.getElementsByTagName('name').item(0);
    if ~isempty(namexml)
        name = char(namexml.getTextContent);
    else
        name = 'Unknown';
    end
    % Get the description of the element
    if ~isempty(node.getElementsByTagName('description').item(0))
        description = char(node.getElementsByTagName('description').item(0).getTextContent);
    end
    % Try to match the style either to the style_has or read the style from
    % the node itself
    style_url = node.getElementsByTagName('styleUrl').item(0);
    if ~isempty(style_url)
        id = char(style_url.getTextContent);
        s = style_hash(id(2:end));
        pointcolor = s.pointcolor; linecolor = s.linecolor; polycolor = s.polycolor;
    else
        [pointcolor,linecolor,polycolor] = parseStyle(node);
    end
    % Find the number of features in this placemark
    number_features = node.getElementsByTagName('coordinates').getLength();
    kmlStructs = cell([1 number_features]);
    count = 1;
    % Handle Points
    points = node.getElementsByTagName('Point');
    for i = 0:points.getLength()-1
        coords = char(points.item(i).getElementsByTagName('coordinates').item(0).getTextContent);
        [Lat,Lon] = parseCoordinates(coords);
        kmlStructs{count}.Geometry = 'Point';
        kmlStructs{count}.Name = name;
        if exist('description','var')
            kmlStructs{count}.Description = description;
        else
            kmlStructs{count}.Description = '';
        end
        kmlStructs{count}.Lon = Lon;
        kmlStructs{count}.Lat = Lat;
        kmlStructs{count}.BoundingBox = [min(Lon) min(Lat);max(Lon) max(Lat)];
        kmlStructs{count}.Folder = folder;
        kmlStructs{count}.Color = pointcolor;
        count = count + 1;
    end
    % Handle Polygons
    polygons = node.getElementsByTagName('Polygon');
    for i = 0:polygons.getLength()-1
        coords = char(polygons.item(i).getElementsByTagName('coordinates').item(0).getTextContent);
        [Lat,Lon] = parseCoordinates(coords);
        kmlStructs{count}.Geometry = 'Polygon';
        kmlStructs{count}.Name = name;
        if exist('description','var')
            kmlStructs{count}.Description = description;
        else
            kmlStructs{count}.Description = '';
        end
        kmlStructs{count}.Lon = [Lon;NaN]';
        kmlStructs{count}.Lat = [Lat;NaN]';
        kmlStructs{count}.BoundingBox = [min(Lon) min(Lat);max(Lon) max(Lat)];
        kmlStructs{count}.Folder = folder;
        kmlStructs{count}.Color = polycolor;
        count = count + 1;
    end
    % Handle Lines
    lines = node.getElementsByTagName('LineString');
    for i = 0:lines.getLength()-1
        coords = char(lines.item(i).getElementsByTagName('coordinates').item(0).getTextContent);
        [Lat,Lon] = parseCoordinates(coords);
        kmlStructs{count}.Geometry = 'Line';
        kmlStructs{count}.Name = name;
        if exist('description','var')
            kmlStructs{count}.Description = description;
        else
            kmlStructs{count}.Description = '';
        end
        kmlStructs{count}.Lon = Lon';
        kmlStructs{count}.Lat = Lat';
        kmlStructs{count}.BoundingBox = [min(Lon) min(Lat);max(Lon) max(Lat)];
        kmlStructs{count}.Folder = folder;
        kmlStructs{count}.Color = linecolor;
        count = count + 1;
    end
    % Compile answers
    kmlStruct = horzcat(kmlStructs{:});
end
function [Lat,Lon] = parseCoordinates(string)
% Read in the "Coordinates" tag of a KML file as quickly as possible.
%
% :param string: Character array from the "Coordinates" tag.
    % Find the coordinates and convert to doubles
    coords = str2double(regexp(string,'[,\s]+','split'));
    coords(isnan(coords)) = [];
    % Determine if these are 3D or 2D coords and shape correctly
    [m,n] = size(coords);
    if length(coords) == sum(string==',') * 2
        coords = reshape(coords,2,m*n/2)';
    else
        coords = reshape(coords,3,m*n/3)';
    end
    % Parse coords to Lat and Lon. If Mapping Toolboxes exists use it to
    % clean the polygon slightly.
    if license('test', 'map_toolbox')
        [Lat, Lon] = poly2ccw(coords(:,2),coords(:,1));
    else
        Lat=coords(:,2);
        Lon=coords(:,1);
    end
end
function [pointcolor,linecolor,polycolor] = parseStyle(node)
% Load style from a "Style" KML tag. If we fail return the color blue.
%
% :param node: A style node to find colors in.
    % Create anon function to get color by style
    get_color_by_style = @(type)char(node.getElementsByTagName(type).item(0).getElementsByTagName('color').item(0).getTextContent);
    
    % Attempt to load color for point geometry
    try
        pointcolorhex = get_color_by_style('IconStyle');
        pointcolor = hex2rgb_kmlwrapper(pointcolorhex);
    catch
        pointcolor = [0.6758    0.8438    0.8984];
    end
    % Attempt to load color for line geometry
    try
        linecolorhex = get_color_by_style('LineStyle');
        linecolor = hex2rgb_kmlwrapper(linecolorhex);
    catch
        linecolor = [0.6758    0.8438    0.8984];
    end
    % Attempt to load color for polygon geometry
    try
        polycolorhex = get_color_by_style('PolyStyle');
        polycolor = hex2rgb_kmlwrapper(polycolorhex);
    catch
        polycolor = [0.6758    0.8438    0.8984];
    end
end
function dec_rgb = hex2rgb_kmlwrapper(hex)
% This code is modified from Chad Greene's original hex2rgb function.
% link: https://www.mathworks.com/matlabcentral/fileexchange/45727-hex2rgb
% Copyright (c) 2014, Chad Greene
% All rights reserved.
%
% Redistribution and use in source and binary forms, with or without
% modification, are permitted provided that the following conditions are met:
%
% * Redistributions of source code must retain the above copyright notice, this
%   list of conditions and the following disclaimer.
% * Redistributions in binary form must reproduce the above copyright notice,
%   this list of conditions and the following disclaimer in the documentation
%   and/or other materials provided with the distribution
% * Neither the name of The University of Texas at Austin nor the names of its
%   contributors may be used to endorse or promote products derived from this
%   software without specific prior written permission.
% THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
% AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
% IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
% DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE
% FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
% DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
% SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
% CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
% OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
% OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
    hex_rgb = hex([7 8 5 6 3 4]);
    dec_rgb = reshape(sscanf(hex_rgb.','%2x'),3,[]).'/255;
end