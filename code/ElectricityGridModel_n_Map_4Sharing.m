%% * Electricity network topology plotting * 
%%%%%%% The figure serves only to illustrate the electricity network topology only %%%%%%%


clear all;
close all;

%%% Read Singapore planning area shape file
disp('Spatial mapping ...');
shpdir = '../data/mp14-plng-area-no-sea-planning-area/MP14_PLNG_AREA_NO_SEA_PL.shp';
disp(shpdir);
sing = shaperead(shpdir);
info = shapeinfo(shpdir);
p1 = info.CoordinateReferenceSystem;
A = readmatrix("../data/PowerGrid/A_Topology_Singapore.csv");
centroids = readmatrix("../data/PowerGrid/centroids_locations.csv");

textString = "Electricity network";  




dataValues = zeros(55, 1);
index_GasGenerator = zeros(1, 3);



for i = 1:length(sing)
    x1 = [sing(i).X];
    y1 = [sing(i).Y];
    [lat,lon] = projinv(p1,x1,y1);
    nanIndices = find(isnan(lon) | isnan(lat));
    startIndex = 1;

    for partIndex = 1:length(nanIndices)
        endIndex = nanIndices(partIndex) - 1;
        
        patch(lon(startIndex:endIndex), lat(startIndex:endIndex), [0.8, 1, 0.8], ...
        'EdgeColor', 'black', ... 
          'LineWidth', 0.5, ... 
          'FaceAlpha', 0.5); 
        
        startIndex = nanIndices(partIndex) + 1;
    end


    if startIndex <= length(lon)
        patch(lon(startIndex:end), lat(startIndex:end), [0.8, 1, 0.8],...
        'EdgeColor', 'black', ... 
          'LineWidth', 0.5, ... 
          'FaceAlpha', 0.5); 
    end
end

hold on;


scatter(centroids(25,1), centroids(25,2), 100, "filled", "hexagram","MarkerEdgeColor","none", "MarkerFaceColor","red");
scatter(centroids(48,1), centroids(48,2), 100, "filled", "hexagram","MarkerEdgeColor","none", "MarkerFaceColor","red");
scatter(centroids(49,1), centroids(49,2), 100, "filled", "hexagram","MarkerEdgeColor","none", "MarkerFaceColor","red");

[n, ~] = size(centroids); 
colors4colormap = flipud(hot); 
cmap = colormap(cool);

line_count = 1;
for i = 1:n
    for j = i+1:n 
    if A(i,j) == 1 

        colorIndex = 103;
        lineColor = cmap(colorIndex, :);
        if ismember(line_count, [74])
            lineColor = [1,0,0];
        end
        quiver(centroids(i,1), centroids(i,2), ...
           centroids(j,1) - centroids(i,1), ...
           centroids(j,2) - centroids(i,2), ...
           0, 'MaxHeadSize', 0.8, 'Color', lineColor,...
           'LineStyle', "-",...
           'LineWidth', 2, 'AutoScale', 'on', 'AutoScaleFactor', 4);
        
          
        mid_x = (centroids(i,1) + centroids(j,1)) / 2;
        mid_y = (centroids(i,2) + centroids(j,2)) / 2;
        

        text(mid_x, mid_y, string(line_count ), ...
         'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', ...
         'Color', 'k'); 
        
        line_count = line_count +1;
    end
    
    end

end

axis off;
  
scatter(103.95, 1.2, 100, "filled", "hexagram","MarkerEdgeColor","none", "MarkerFaceColor","red");
text(103.90, 1.175, "Fossil fuel generator", "FontSize", 12);


%{  %}

%% ---------------- Export (vector PDF) ----------------
outPath = "../results/TopologySingaporeGridVisualisation_noArrow_" +...
string(datetime("now", "Format", "yyyy-MM-dd")) + ".pdf"

exportgraphics(gcf, outPath, ...
'ContentType','vector', 'BackgroundColor','none');

%{ %}

