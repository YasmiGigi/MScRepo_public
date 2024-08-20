%%
clear all, clc, close all, home

%% load 3D EBSD data

fileName        = 'data\IN718new.dream3d';
ebsd            = import_dream3D(fileName);

%% rotational

deg             = 0:180;
vector          = [xvector,yvector,zvector];
len             = length(deg)*length(vector);

diagonal        = ceil(sqrt(max(ebsd.x)^2+max(ebsd.y)^2+max(ebsd.z)^2));
[x, y, z]       = meshgrid(-diagonal/2:diagonal/2,-diagonal/2:diagonal/2,0);
plane           = vector3d(x(:),y(:),z(:));

%% variables

% data filtering and correction
oriRemove       = orientation.byEuler(270*degree,0,0);
iqRemove        = 8000;
grainRemove     = 15;
voxelSize       = 1.5;

% grain boundaries
lowAngle        = 2*degree;
highAngle       = 15*degree;
smoothIt        = 20;

% kernel halfwidth
delta           = 7*degree;

% size binning thresholds
thresholdSmall  = 300;
thresholdLarge  = 1800;

% size "phases"
small           = crystalSymmetry('m-3m',...
                    'mineral', 'small',  'color', [255 161 170]./255);
medium          = crystalSymmetry('m-3m',...
                    'mineral', 'medium', 'color', [190 238 175]./255);
large           = crystalSymmetry('m-3m',...
                    'mineral', 'large',  'color', [180 203 240]./255);
sizes           = {small, medium, large};
for s = 3:5
    ebsd.CSList{s}   = sizes{s-2};
    ebsd.phaseMap(s) = s-1;
end

% volume
odf_3D          = calcDensity(ebsd.orientations,'halfwidth',7*degree);
[~,ori]         = max(odf_3D);
cube            = orientation.byEuler(ori.phi1,ori.Phi,ori.phi2,ebsd.CS);
f               = fibre(Miller(0,0,1,ebsd.CS),zvector);
volRadius       = 20*degree;

% counting progress
k               = 1;

%% variables for saving results

%initialize variables
rotationSave    = [];                               % rotation
mapArea         = [];                               % area size of EBSD map
grainSizeRatio  = [];                               % ratio of grain sizes
noGrains        = zeros([len,4,1]);                 % number of grains
shapeFactorStat = zeros([len,3,4]);                 % shape factor statitics
areaStat        = zeros([len,3,4]);                 % grain area statistics
textureIndex    = zeros([len,4,1]);                 % texture index
volumeCube      = zeros([len,4,1]);                 % volume, texture
volumeFiber     = zeros([len,4,1]);                 % volume, texture
boundDen        = zeros([len,3,4]);                 % boundary density
GOS             = zeros([len,3,4]);                 % grain orientation spread


% subgrains
noGrainsSub     = zeros([len,4,1]); 
SFStatSub       = zeros([len,3,4]); 
areaStatSub     = zeros([len,3,4]);

%% orientation analysis

for i = 1:length(vector)

    % rotational vector
    v = vector(:,i);

    if      i == 1
        rotZ = rotation.byAxisAngle(xvector,90*degree);
    elseif  i == 2
        rotZ = rotation.byAxisAngle(yvector,90*degree);
    elseif  i == 3
        rotZ = rotation.byAxisAngle(yvector,90*degree);
    end

    for j = 1:length(deg)

        % print loop progress
        count = append('progress: ',num2str(k),'/',num2str(length(deg)*3));
        clc;disp('rotational');disp(count);k=k+1;

        % rotation of selected slice, save rotation for future reference
        rot = rotation.byAxisAngle(v,deg(j)*degree);
        rotationSave = [rotationSave; [rot.phi1 rot.Phi rot.phi2]/degree];

        % reference frame to rotate for grain calculation
        rotReference1 = rotation.byEuler(-rot.phi2,-rot.Phi,-rot.phi1);
        rotReference2 = rotation.byEuler(-rotZ.phi2,-rotZ.Phi,-rotZ.phi1);

        % rotate empty grid, move to center of cuboid
        newGrid     = rotate(plane,rotZ);
        newGrid     = rotate(newGrid,rot);
        newGrid.x   = newGrid.x+max(ebsd.x)/2;
        newGrid.y   = newGrid.y+max(ebsd.y)/2;
        newGrid.z   = newGrid.z+max(ebsd.z)/2;

        % interpolate data from 3D EBSD data based on selected slice
        ebsdSlice   = ebsd.interp(newGrid);

        % rotate reference frame, keep Euler angles
        ebsdSlice   = rotate(ebsdSlice,rotReference1,'keepEuler');
        ebsdSlice   = rotate(ebsdSlice,rotReference2,'keepEuler');

        % correction of map size
        ebsdSlice.pos = ebsdSlice.pos*voxelSize;

        % filter grain data
        ebsdSlice(ebsdSlice.orientations==oriRemove) = 'notIndexed';
        ebsdSlice(ebsdSlice.iq<=iqRemove)            = 'notIndexed';

        % convert EBSD3 to EBSD to allow calcGrains()
        ebsdSlice = EBSD(ebsdSlice);

        % calculate grains, remove small grains
        [grains, ebsdSlice.grainId, ebsdSlice.mis2mean] = calcGrains(ebsdSlice,'angle',highAngle);
        ebsdSlice(grains(grains.grainSize <= grainRemove)) = 'notIndexed';

        % calculate grains and subgrains, smooth grain boundaries
        [grains, ebsdSlice.grainId, ebsdSlice.mis2mean] = calcGrains(ebsdSlice,'angle',[lowAngle, highAngle]);
        grains = smooth(grains,smoothIt,'moveTriplePoints');

        % remove non-indexed EBSD data and grains
        ebsdSlice   = ebsdSlice('indexed');
        grains      = grains('indexed');

        % save map size
        xSize       = (max(ebsdSlice.x)-min(ebsdSlice.x))./voxelSize;
        ySize       = (max(ebsdSlice.y)-min(ebsdSlice.y))./voxelSize;
        mapArea     = [mapArea; xSize*ySize xSize ySize];

        % grain binned based on sizes
        grains.phase(grains.area<=thresholdSmall) = 2;
        grains.phase(grains.area> thresholdSmall) = 3;
        grains.phase(grains.area> thresholdLarge) = 4;

        % grain size ratios
        grainSizeRatio   = [grainSizeRatio; sum(grains('small').area)/sum(grains.area) ...
                                            sum(grains('medium').area)/sum(grains.area) ...
                                            sum(grains('large').area)/sum(grains.area)];

        % find no. of grains
        noGrains(j+(length(deg)*(i-1)),:)   = [length(grains) ...               % all grains
                            length(grains('small')) ...     % small grains
                            length(grains('medium')) ...    % medium grains
                            length(grains('large'))];       % large grains

        % find shape factor
        shapeFactorCalc         = shapeFactor(grains);
        shapeFactorStat(j+(length(deg)*(i-1)),:,1)  = [mean(shapeFactorCalc) ...        % all grains
                                    median(shapeFactorCalc) ...
                                    std(shapeFactorCalc)];
        shapeFactorCalc         = shapeFactor(grains('small'));
        shapeFactorStat(j+(length(deg)*(i-1)),:,2)  = [mean(shapeFactorCalc) ...        % small grains      
                                    median(shapeFactorCalc) ...
                                    std(shapeFactorCalc)];
        shapeFactorCalc         = shapeFactor(grains('medium'));
        shapeFactorStat(j+(length(deg)*(i-1)),:,3)  = [mean(shapeFactorCalc) ...        % medium grains
                                    median(shapeFactorCalc) ...
                                    std(shapeFactorCalc)];
        shapeFactorCalc         = shapeFactor(grains('large'));
        shapeFactorStat(j+(length(deg)*(i-1)),:,4)  = [mean(shapeFactorCalc) ...        % large grains
                                median(shapeFactorCalc) ...
                                std(shapeFactorCalc)];

        % find grain area stats
        areaCalc                = grains.area;
        areaStat(j+(length(deg)*(i-1)),:,1)         = [mean(areaCalc) ...               % all grains
                                    median(areaCalc) ...
                                    std(areaCalc)];
        areaCalc                = grains('small').area;
        areaStat(j+(length(deg)*(i-1)),:,2)         = [mean(areaCalc) ...               % small grains
                                    median(areaCalc) ...
                                    std(areaCalc)];
        areaCalc                = grains('medium').area;
        areaStat(j+(length(deg)*(i-1)),:,3)         = [mean(areaCalc) ...               % medium grains
                                    median(areaCalc) ...
                                    std(areaCalc)];
        areaCalc                = grains('large').area;
        areaStat(j+(length(deg)*(i-1)),:,4)         = [mean(areaCalc) ...               % large grains
                                    median(areaCalc) ...
                                    std(areaCalc)];
        % find texture index
        odf                     = calcDensity(ebsdSlice.orientations,'halfwidth',delta);
        odfSmall                = calcDensity(ebsdSlice(grains('small')).orientations,'halfwidth',delta);
        odfMedium               = calcDensity(ebsdSlice(grains('medium')).orientations,'halfwidth',delta);
        odfLarge                = calcDensity(ebsdSlice(grains('large')).orientations,'halfwidth',delta);
    
        textureIndex(j+(length(deg)*(i-1)),:) = [norm(odf)^2 norm(odfSmall)^2 norm(odfMedium)^2 norm(odfLarge)^2];
    
        % volume
        vol                     = volume(odf,cube,volRadius);
        volSmall                = volume(odfSmall,cube,volRadius);
        volMedium               = volume(odfMedium,cube,volRadius);
        volLarge                = volume(odfLarge,cube,volRadius);
    
        volumeCube(j+(length(deg)*(i-1)),:)  = [vol volSmall volMedium volLarge];
    
        vol                     = volume(odf,f,volRadius);
        volSmall                = volume(odfSmall,f,volRadius);
        volMedium               = volume(odfMedium,f,volRadius);
        volLarge                = volume(odfLarge,f,volRadius);

        volumeFiber(j+(length(deg)*(i-1)),:) = [vol volSmall volMedium volLarge];


        % subgrain boundary density
    
        boundDen(j+(length(deg)*(i-1)),:,1) = [mean(grains.subBoundaryLength./grains.area) ...
                                        median(grains.subBoundaryLength./grains.area) ...
                                        std(grains.subBoundaryLength./grains.area)];
        boundDen(j+(length(deg)*(i-1)),:,2) = [mean(grains('small').subBoundaryLength./grains('small').area) ...
                                        median(grains('small').subBoundaryLength./grains('small').area) ...
                                        std(grains('small').subBoundaryLength./grains('small').area)];
        boundDen(j+(length(deg)*(i-1)),:,3) = [mean(grains('medium').subBoundaryLength./grains('medium').area) ...
                                        median(grains('medium').subBoundaryLength./grains('medium').area) ...
                                        std(grains('medium').subBoundaryLength./grains('medium').area)];
        boundDen(j+(length(deg)*(i-1)),:,4) = [mean(grains('large').subBoundaryLength./grains('large').area) ...
                                        median(grains('large').subBoundaryLength./grains('large').area) ...
                                        std(grains('large').subBoundaryLength./grains('large').area)];
    
        % grain orientation spread
        GOS(j+(length(deg)*(i-1)),:,1)      = [mean(grains.GOS)./degree ...                 % all grains
                                        median(grains.GOS)./degree ...
                                        std(grains.GOS)./degree];
        GOS(j+(length(deg)*(i-1)),:,2)      = [mean(grains('small').GOS)./degree ...        % small grains
                                        median(grains('small').GOS)./degree ...
                                        std(grains('small').GOS)./degree];
        GOS(j+(length(deg)*(i-1)),:,3)      = [mean(grains('medium').GOS)./degree ...       % medium grains
                                        median(grains('medium').GOS)./degree ...
                                        std(grains('medium').GOS)./degree];
        GOS(j+(length(deg)*(i-1)),:,4)      = [mean(grains('large').GOS)./degree ...        % large grains
                                        median(grains('large').GOS)./degree ...
                                        std(grains('large').GOS)./degree];
        
        % subgrain boundary length
        subBoundLength(j+(length(deg)*(i-1)),:)     = [sum(grains.subBoundaryLength) ...             % all grains
                                    sum(grains('small').subBoundaryLength) ...   % small grains
                                    sum(grains('medium').subBoundaryLength) ...  % medium grains
                                    sum(grains('large').subBoundaryLength)];    % large grains
    
        ebsdSlice(grains('small'))  = 'notIndexed';
        ebsdSlice(grains('medium')) = 'notIndexed';
        
        % subgrains
        [subgrains, ebsdSlice.grainId] = calcGrains(ebsdSlice,'angle',lowAngle);
        subgrains = smooth(subgrains,smoothIt);
        subgrains = subgrains('indexed');
    
        % subgrains binned based on sizes
        subgrains.phase(subgrains.area<=thresholdSmall) = 2;
        subgrains.phase(subgrains.area> thresholdSmall) = 3;
        subgrains.phase(subgrains.area> thresholdLarge) = 4;
        
        % find no. of subgrains
        noGrainsSub(j+(length(deg)*(i-1)),:)        = [length(subgrains) ...            % all grains
                                    length(subgrains('small')) ...  % small grains
                                    length(subgrains('medium')) ... % medium grains
                                    length(subgrains('large'))];    % large grains
    
        % find shape factor
        shapeFactorCalc         = shapeFactor(subgrains);
        SFStatSub(j+(length(deg)*(i-1)),:,1)        = [mean(shapeFactorCalc) ...        % all grains
                                    median(shapeFactorCalc) ...
                                    std(shapeFactorCalc)];
        shapeFactorCalc         = shapeFactor(subgrains('small'));
        SFStatSub(j+(length(deg)*(i-1)),:,2)        = [mean(shapeFactorCalc) ...        % small grains      
                                    median(shapeFactorCalc) ...
                                    std(shapeFactorCalc)];
        shapeFactorCalc         = shapeFactor(subgrains('medium'));
        SFStatSub(j+(length(deg)*(i-1)),:,3)        = [mean(shapeFactorCalc) ...        % medium grains
                                    median(shapeFactorCalc) ...
                                    std(shapeFactorCalc)];
        shapeFactorCalc         = shapeFactor(subgrains('large'));
        SFStatSub(j+(length(deg)*(i-1)),:,4)        = [mean(shapeFactorCalc) ...        % large grains
                                    median(shapeFactorCalc) ...
                                std(shapeFactorCalc)];

        % find subgrain area stats
        areaCalc                = subgrains.area;
        areaStatSub(j+(length(deg)*(i-1)),:,1)      = [mean(areaCalc) ...               % all grains
                                    median(areaCalc) ...
                                    std(areaCalc)];
        areaCalc                = subgrains('small').area;
        areaStatSub(j+(length(deg)*(i-1)),:,2)      = [mean(areaCalc) ...               % small grains
                                    median(areaCalc) ...
                                    std(areaCalc)];
        areaCalc                = subgrains('medium').area;
        areaStatSub(j+(length(deg)*(i-1)),:,3)      = [mean(areaCalc) ...               % medium grains
                                    median(areaCalc) ...
                                    std(areaCalc)];
        areaCalc                = subgrains('large').area;
        areaStatSub(j+(length(deg)*(i-1)),:,4)      = [mean(areaCalc) ...               % large grains
                                    median(areaCalc) ...
                                    std(areaCalc)];

    end
end

%% save results
results = array2table([rotationSave mapArea grainSizeRatio...
    noGrains ...
    shapeFactorStat(:,:,1) shapeFactorStat(:,:,2) shapeFactorStat(:,:,3) shapeFactorStat(:,:,4) ...
    areaStat(:,:,1) areaStat(:,:,2) areaStat(:,:,3) areaStat(:,:,4) ...
    textureIndex ...
    volumeCube ...
    volumeFiber ...
    boundDen(:,:,1) boundDen(:,:,2) boundDen(:,:,3) boundDen(:,:,4) ...
    GOS(:,:,1) GOS(:,:,2) GOS(:,:,3) GOS(:,:,4)],'VariableNames',{'phi1', 'Phi','phi2',...
    'Map Area',                     'x size',                       'y size',...
    'Ratio of Small Grains',        'Ratio of Medium Grains',       'Ratio of Large Grains',...
    'No. of Grains',                'No. of Grains (Small)',        'No. of Grains (Medium)',                   'No. of Grains (Large)',...
    'Mean Shape Factor',            'Median Shape Factor',          'Standard Deviation Shape Factor',...
    'Mean Shape Factor (Small)',    'Median Shape Factor (Small)',  'Standard Deviation Shape Factor (Small)',...
    'Mean Shape Factor (Medium)',   'Median Shape Factor (Medium)', 'Standard Deviation Shape Factor (Medium)',...
    'Mean Shape Factor (Large)',    'Median Shape Factor (Large)',  'Standard Deviation Shape Factor (Large)',...
    'Mean Grain Area',              'Median Grain Area',            'Standard Deviation Grain Area',...
    'Mean Grain Area (Small)',      'Median Grain Area (Small)',    'Standard Deviation Grain Area (Small)',...
    'Mean Grain Area (Medium)',     'Median Grain Area (Medium)',   'Standard Deviation Grain Area (Medium)',...
    'Mean Grain Area (Large)',      'Median Grain Area (Large)',    'Standard Deviation Grain Area (Large)',...
    'Texture Index',                'Texture Index (Small)',        'Texture Index (Medium)',                   'Texture Index (Large)',...
    'Volume, Cube',                 'Volume, Cube (Small)',         'Volume, Cube (Medium)',                    'Volume, Cube (Large)',...
    'Volume, Fiber',                'Volume, Fiber (Small)',        'Volume, Fiber (Medium)',                   'Volume, Fiber (Large)',...
    'Mean Subgrain Boundary Density',           'Median Boundary Density',                  'Standard Deviation Boundary Density',...
    'Mean Subgrain Boundary Density (Small)',   'Median Boundary Density (Small)',          'Standard Deviation Boundary Density (Small)',...
    'Mean Subgrain Boundary Density (Medium)',  'Median Boundary Density (Medium)',         'Standard Deviation Boundary Density (Medium)',...
    'Mean Subgrain Boundary Density (Large)',   'Median Boundary Density (Large)',          'Standard Deviation Boundary Density (Large)',...
    'Mean GOS',                     'Median GOS',                   'Standard Deviation GOS',...
    'Mean GOS (Small)',             'Median GOS (Small)',           'Standard Deviation GOS (Small)',...
    'Mean GOS (Medium)',            'Median GOS (Medium)',          'Standard Deviation GOS (Medium)',...
    'Mean GOS (Large)',             'Median GOS (Large)',           'Standard Deviation GOS (Large)'});

% save all results in .csv file
writetable(results,'output\results_IN718_rotational.csv')

disp('results saved')

%% save subgrain results

resultsSub = array2table([rotationSave ...
    noGrainsSub ...
    boundDen(:,:,1) boundDen(:,:,2) boundDen(:,:,3) boundDen(:,:,4) ...
    GOS(:,:,1) GOS(:,:,2) GOS(:,:,3) GOS(:,:,4) ...
    subBoundLength ...
    SFStatSub(:,:,1) SFStatSub(:,:,2) SFStatSub(:,:,3) SFStatSub(:,:,4) ...
    areaStatSub(:,:,1) areaStatSub(:,:,2) areaStatSub(:,:,3) areaStatSub(:,:,4)],'VariableNames',{'phi1', 'Phi','phi2',...
    'No. of Grains (Subgrains)',                'No. of Grains (Subgrains, Small)',         'No. of Grains (Subgrains, Medium)',        'No. of Grains (Subgrains, Large)',...
    'Mean Subgrain Boundary Density',           'Median Boundary Density',                  'Standard Deviation Boundary Density',...
    'Mean Subgrain Boundary Density (Small)',   'Median Boundary Density (Small)',          'Standard Deviation Boundary Density (Small)',...
    'Mean Subgrain Boundary Density (Medium)',  'Median Boundary Density (Medium)',         'Standard Deviation Boundary Density (Medium)',...
    'Mean Subgrain Boundary Density (Large)',   'Median Boundary Density (Large)',          'Standard Deviation Boundary Density (Large)',...
    'Mean GOS',                                 'Median GOS',                               'Standard Deviation GOS',...
    'Mean GOS (Small)',                         'Median GOS (Small)',                       'Standard Deviation GOS (Small)',...
    'Mean GOS (Medium)',                        'Median GOS (Medium)',                      'Standard Deviation GOS (Medium)',...
    'Mean GOS (Large)',                         'Median GOS (Large)',                       'Standard Deviation GOS (Large)',...
    'Subgrain Boundary Length',                 'Subgrain Boundary Length (Small)',         'Subgrain Boundary Length (Medium)',         'Subgrain Boundary Length (Large)',...
    'Mean Shape Factor (Subgrains)',            'Median Shape Factor (Subgrains)',          'Standard Deviation Shape Factor (Subgrains)',...
    'Mean Shape Factor (Subgrains, Small)',     'Median Shape Factor (Subgrains, Small)',   'Standard Deviation Shape Factor (Subgrains, Small)',...
    'Mean Shape Factor (Subgrains, Medium)',    'Median Shape Factor (Subgrains, Medium)',  'Standard Deviation Shape Factor (Subgrains, Medium)',...
    'Mean Shape Factor (Subgrains, Large)',     'Median Shape Factor (Subgrains, Large)',   'Standard Deviation Shape Factor (Subgrains, Large)',...
    'Mean Grain Area (Subgrains)',              'Median Grain Area (Subgrains)',            'Standard Deviation Grain Area (Subgrains)',...
    'Mean Grain Area (Subgrains, Small)',       'Median Grain Area (Subgrains, Small)',     'Standard Deviation Grain Area (Subgrains, Small)',...
    'Mean Grain Area (Subgrains, Medium)',      'Median Grain Area (Subgrains, Medium)',    'Standard Deviation Grain Area (Subgrains, Medium)',...
    'Mean Grain Area (Subgrains, Large)',       'Median Grain Area (Subgrains, Large)',     'Standard Deviation Grain Area (Subgrains, Large)'});

% save all results in .csv file
writetable(resultsSub,'output\results_IN718_rotational_subgrains.csv')

disp('results (subgrains) saved')