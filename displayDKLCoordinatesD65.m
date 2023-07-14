% ============================================================================
% function displayDKLCoordinatesD65 
% SourceCode: displayDKLCoordinates.m 
% [rgb0,rgb1,nRGB0,nRGB1] = displayDKLCoordinatesD65(kdlTheta)
% plots colors defined in the DKL space in standard color spaces
% and returns RGB values of the corresponding DKL stimulus
% rgb0 - computed using Lablib way
% rgb1 - computed using Matlab (XYZ to RGB conversion)
% nRGB0 and nRGB1 - normalized rgb0 and rgb1
% ============================================================================

function [rgb0,rgb1,nRGB0,nRGB1] = displayDKLCoordinatesD65(kdlTheta)

if ~exist('kdlTheta','var');            kdlTheta=0;                     end

% x and y coordinates of the primaries, computed from the spectra measured using PR-655
% Colors were displayed using Lablib in a gamma corrected monitor (BenQXL2411) 
% after putting the x and y coordinates of the primaries.
% These coordinates x and y coordinates were computed from the XYZ values 
% of the ICC profile made using i1DisplayPro. 

CIEx.r = 0.6432; 
CIEx.g = 0.3112; 
CIEx.b = 0.1522;

CIEy.r = 0.3299;
CIEy.g = 0.6104;
CIEy.b = 0.0584;

% Standard CIE coordionates of the white points:
CIEx.wEE = 1/3; % Equal Energy
CIEy.wEE = 1/3;

CIEx.wD65 = 0.31272;  % D65 white point 
CIEy.wD65 = 0.32902;


% Enter the x and y coordinates of the stimulus
readParentPath = fullfile(pwd,'Measurements','Rig1Display');
fileNames      = 'Rig1DisplayProfilingPR655_cieCoordinates';
cieStimsAll = readtable(fullfile(readParentPath,fileNames));

%rgbWBG
startPosPri=76; endPosPri=81;
priWhite_x =  table2array(cieStimsAll(8,startPosPri:endPosPri));
priWhite_y =  table2array(cieStimsAll(9,startPosPri:endPosPri));
priWhite_Y =  table2array(cieStimsAll(3,startPosPri:endPosPri));
priWhite_Y_norm = priWhite_Y/max(priWhite_Y);

% 16 DKL stimulus now, azimuth 0:22.5:33.7.5
startPosDKL=46; endPosDKL=61;
stimCIEx = table2array(cieStimsAll(8,startPosDKL:endPosDKL));
stimCIEy = table2array(cieStimsAll(9,startPosDKL:endPosDKL));
stimCIEY = table2array(cieStimsAll(3,startPosDKL:endPosDKL));
normStimCIEY = stimCIEY/max(priWhite_Y);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% xy Coordinates
% subplot(231);
subplot(231);
axis square;
cieCoordinates=load('cieCoordinates.mat');
plot(cieCoordinates.cieCoordinates(:,2), cieCoordinates.cieCoordinates(:,3),'k'); axis([0 1 0 1]);
hold on;

% plotting 
plot(CIEx.r,CIEy.r,'marker','o','color','r');
plot(CIEx.g,CIEy.g,'marker','o','color','g');
plot(CIEx.b,CIEy.b,'marker','o','color','b');
plot(CIEx.wEE,CIEy.wEE,'marker','o','color',[0.5 0.5 0.5]);

%plot the observed primaries and white point:

plot(priWhite_x(1),priWhite_y(1),'marker','*','color','r','markerfacecolor','r');
plot(priWhite_x(2),priWhite_y(2),'marker','*','color','g','markerfacecolor','g');
plot(priWhite_x(3),priWhite_y(3),'marker','*','color','b','markerfacecolor','b');
plot(priWhite_x(4),priWhite_y(4),'marker','*','color',[0.5 0.5 0.5],'markerfacecolor',[0.5 0.5 0.5]); % plotted for white as of now

% 
line([CIEx.b CIEx.g],[CIEy.b CIEy.g],'color','k');
line([CIEx.b CIEx.r],[CIEy.b CIEy.r],'color','k');
line([CIEx.g CIEx.r],[CIEy.g CIEy.r],'color','k');
xlabel('x'); ylabel('y');
xlim([0 1]); ylim([0 1]);
title('CIE xy');

% xyY Coordinates
% subplot(234);
subplot(233);
% Interestingly, once you define the white point, the fraction of luminance
% contributed by the 3 phosphors gets fixed. Note that the second row of
% this Matrix gives the equalEnergy colors in Lablib
Mee  = RGBToXYZMatrix(CIEx.r, CIEy.r, CIEx.g, CIEy.g, CIEx.b, CIEy.b, CIEx.wEE, CIEx.wEE);
Md65 = RGBToXYZMatrix(CIEx.r, CIEy.r, CIEx.g, CIEy.g, CIEx.b, CIEy.b, CIEx.wD65, CIEy.wD65);

% EE 
XYZrEE = Mee*[1 0 0]'; xyYrEE = XYZToxyY(XYZrEE);
XYZgEE = Mee*[0 1 0]'; xyYgEE = XYZToxyY(XYZgEE);
XYZbEE = Mee*[0 0 1]'; xyYbEE = XYZToxyY(XYZbEE);
XYZwEE = Mee*[0.5 0.5 0.5]'; xyYwEE = XYZToxyY(XYZwEE);

% DD 
XYZrDD = Md65*[1 0 0]'; xyYrDD = XYZToxyY(XYZrDD);
XYZgDD = Md65*[0 1 0]'; xyYgDD = XYZToxyY(XYZgDD);
XYZbDD = Md65*[0 0 1]'; xyYbDD = XYZToxyY(XYZbDD);
XYZbDD = Md65*[0.5 0.5 0.5]'; xyYwDD = XYZToxyY(XYZbDD );

% predicted primaries and white point:
plot3(xyYrDD(1),xyYrDD(2),xyYrDD(3)*max(priWhite_Y),'marker','o','color','r'); hold on;
plot3(xyYgDD(1),xyYgDD(2),xyYgDD(3)*max(priWhite_Y),'marker','o','color','g');
plot3(xyYbDD(1),xyYbDD(2),xyYbDD(3)*max(priWhite_Y),'marker','o','color','b');
plot3(xyYwDD(1),xyYwDD(2),xyYwDD(3)*max(priWhite_Y),'marker','o','color',[0.5 0.5 0.5]);

% observed primaries and white point:
plot3(priWhite_x(1),priWhite_y(1),priWhite_Y_norm(1)*max(priWhite_Y),'marker','*','color','r','markerfacecolor','r'); hold on;
plot3(priWhite_x(2),priWhite_y(2),priWhite_Y_norm(2)*max(priWhite_Y),'marker','*','color','g','markerfacecolor','g');
plot3(priWhite_x(3),priWhite_y(3),priWhite_Y_norm(3)*max(priWhite_Y),'marker','*','color','b','markerfacecolor','b');
plot3(priWhite_x(4),priWhite_y(4),priWhite_Y_norm(6)*max(priWhite_Y),'marker','*','color',[0.5 0.5 0.5],'markerfacecolor',[0.5 0.5 0.5]); %plotting for gray

% observed primaries and white point: 2nd set after correction

% connect the dots
line([xyYbDD(1) xyYgDD(1)],[xyYbDD(2) xyYgDD(2)],[xyYbDD(3)*max(priWhite_Y) xyYgDD(3)*max(priWhite_Y)],'color','k');
line([xyYrDD(1) xyYgDD(1)],[xyYrDD(2) xyYgDD(2)],[xyYrDD(3)*max(priWhite_Y) xyYgDD(3)*max(priWhite_Y)],'color','k');
line([xyYrDD(1) xyYbDD(1)],[xyYrDD(2) xyYbDD(2)],[xyYrDD(3)*max(priWhite_Y) xyYbDD(3)*max(priWhite_Y)],'color','k');
xlabel('x'); ylabel('y'); zlabel('Y');
xlim([0 1]); ylim([0 1]); 
zlim([0 120]);
title('CIE xyY - Lablib Colors');

%%%%%%%%%%%%%%%%%%%%%%% Macleod Boynton (MB) Space %%%%%%%%%%%%%%%%%%%%%%%%
subplot(232);
axis square;

% We use the convention used in Lablib
% predicted primaries and white point
[bmr0,bmb0] = xy2MB(CIEx.r,CIEy.r);
[bmr1,bmb1] = xy2MB(CIEx.g,CIEy.g);
[bmr2,bmb2] = xy2MB(CIEx.b,CIEy.b);
[bmr3,bmb3] = xy2MB(CIEx.wEE,CIEy.wEE);

% observed primaries and white point:
[bmr4,bmb4] = xy2MB(priWhite_x(1),priWhite_y(1));
[bmr5,bmb5] = xy2MB(priWhite_x(2),priWhite_y(2));
[bmr6,bmb6] = xy2MB(priWhite_x(3),priWhite_y(3));
[bmr7,bmb7] = xy2MB(priWhite_x(4),priWhite_y(4));


mbLocus = load('MBlocus.mat');
plot(mbLocus.rMB,mbLocus.bMB,'k');
hold on;
% plot(bmr0,bmb0,'marker','o','color','r','markerfacecolor','r');
% plot(bmr1,bmb1,'marker','o','color','g','markerfacecolor','g');
% plot(bmr2,bmb2,'marker','o','color','b','markerfacecolor','b');
% plot(bmr3,bmb3,'marker','o','color',[0.5 0.5 0.5],'markerfacecolor',[0.5 0.5 0.5]);

% plot predicted primaries and white point
plot(bmr0,bmb0,'marker','o','color','r');
plot(bmr1,bmb1,'marker','o','color','g');
plot(bmr2,bmb2,'marker','o','color','b');
plot(bmr3,bmb3,'marker','o','color',[0.5 0.5 0.5]);

% plot observed primaries and white point
plot(bmr4,bmb4,'marker','*','color','r','markerfacecolor','r');
plot(bmr5,bmb5,'marker','*','color','g','markerfacecolor','g');
plot(bmr6,bmb6,'marker','*','color','b','markerfacecolor','b');
plot(bmr7,bmb7,'marker','*','color',[0.5 0.5 0.5],'markerfacecolor',[0.5 0.5 0.5]);

% connect the dots
line([bmr0 bmr1],[bmb0 bmb1],'color','k');
line([bmr0 bmr2],[bmb0 bmb2],'color','k');
line([bmr1 bmr2],[bmb1 bmb2],'color','k');
ylim([0 0.25]);
title('Macleod Boynton');
xlabel('r'); ylabel('b');

%%%%%%%%%%% Compute Carginal Green and Yellow and show them %%%%%%%%%%%%%%%
bmr = [bmr0 bmr1 bmr2 bmr3];
bmb = [bmb0 bmb1 bmb2 bmb3];
calibratedColor = computeKDLColors(bmr,bmb);

% Cardinal Green
bmrCarGreen = calibratedColor.cardinalGreen.green*bmr1+calibratedColor.cardinalGreen.blue*bmr2;
scalingFactor_cb = bmr3 - bmrCarGreen;

% Plot Cardinal Green axis on both MB and xy planes
r1 = bmr3 - scalingFactor_cb;
r2 = bmr3 + scalingFactor_cb;
b = bmb3;
[x1,y1] = MB2xy(r1,b);
[x2,y2] = MB2xy(r2,b);

subplot(231);
line([x1 x2],[y1 y2],'color','k');

subplot(232);
line([r1 r2],[b b],'color','k');

% Cardinal Yellow
bmbCarYellow = calibratedColor.cardinalYellow.red*bmb0+calibratedColor.cardinalYellow.green*bmb1;
scalingFactor_tc = bmb3 - bmbCarYellow;

% Plot the axis
b1 = bmb3 - scalingFactor_tc;
b2 = bmb3 + scalingFactor_tc;
r = bmr3;
[x1,y1] = MB2xy(r,b1);
[x2,y2] = MB2xy(r,b2);

subplot(231);
line([x1 x2],[y1 y2],'color','k');

subplot(232);
line([r r],[b1 b2],'color','k');

%%%%%%%%%%% Find appropriate colors using two approaches %%%%%%%%%%%%%%%%%%

kdlConstants = calibratedColorTokdlConstants(calibratedColor); % Lablib

rgb0 = zeros(16,3); % Colors using Lablib
rgb1 = zeros(16,3); % Colors by directly converting XYZ to RGB

kdlPhiVals = 0:22.5:337.5;
% kdlPhiVals = 0:90:270;

for i=1:length(kdlPhiVals)
%     kdlPhi = (i-1)*10;
    kdlPhi = kdlPhiVals(i);
    % Lablib
    [rgb,lum,cb,tc] = kdlToRGB(kdlConstants,kdlPhi,kdlTheta);
    rgb0(i,:) = [rgb.red rgb.green rgb.blue]; % Actual rgb coordinates using Lablib
    nRGB0 = normalizeColors(rgb,kdlPhi,kdlTheta); % Normalize if needed to keep within range
    
    % Matlab
    % Macleod Boynton coordinates
    bmr = bmr3 + cb*scalingFactor_cb;
    bmb = bmb3 + tc*scalingFactor_tc;
    
    [x,y] = MB2xy(bmr,bmb); % Convert MB to xy
    xyYtmp = [x y 0.5+lum/2]'; % create a color with same xy (and hence MB coordinates) but different luminance
    XYZtmp = xyYToXYZ(xyYtmp);
    rgb1(i,:) = Mee\XYZtmp;
    
    % Show coordinates in xy and MB spaces. Use Lablib colors
    subplot(231); % CIE xy space
    plot(x,y,'marker', 'o','color',[nRGB0.red nRGB0.green nRGB0.blue]); %predicted
    plot(stimCIEx(i),stimCIEy(i),'marker', '*','color',[nRGB0.red nRGB0.green nRGB0.blue]); %predicted

    % actual MB
    [bmrStim,bmbStim] = xy2MB(stimCIEx(i),stimCIEy(i));
    
    subplot(232); % MB space
    plot(bmr,bmb,'marker', 'o','color',[nRGB0.red nRGB0.green nRGB0.blue]); %predicted
    plot(bmrStim,bmbStim,'marker', '*','color',[nRGB0.red nRGB0.green nRGB0.blue]); %actual
%     text(0.8,0.2,'o - Predicted','FontSize',8);
%     text(0.8,0.17,'* - Observed','FontSize',8);
    
    % Show stimuli in normalized DKL space (cb, tc, lum)
    subplot(233)
    plot3(cb,tc,lum,'marker', 'o','color',[nRGB0.red nRGB0.green nRGB0.blue],'markerfacecolor',[nRGB0.red nRGB0.green nRGB0.blue]);
    hold on;
    xlabel('Constant Blue (cb)'); ylabel('Constant R&G (tc)'); zlabel('Luminance (lum)');
    xlim([-1 1]); ylim([-1 1]); zlim([-1 1]);
    
    % Show stimuli in xyY space. Use Lablib colors
    subplot(233);
%     subplot(223);
%     plot3(x,y,0.5+lum/2,'marker', 'o','color',[nRGB0.red nRGB0.green nRGB0.blue],'markerfacecolor',[nRGB0.red nRGB0.green nRGB0.blue]);
    plot3(x,y,(0.5+lum/2)*max(priWhite_Y),'marker', 'o','color',[nRGB0.red nRGB0.green nRGB0.blue]); %predicted
    
    plot3(stimCIEx(i),stimCIEy(i),normStimCIEY(i)*max(priWhite_Y),'marker', '*','color',[nRGB0.red nRGB0.green nRGB0.blue]); %actual
     
    subplot(234);
    plot(kdlPhiVals(i),(0.5+lum/2)*max(priWhite_Y),'marker', 'o','color',[nRGB0.red nRGB0.green nRGB0.blue]);   hold on
    plot(kdlPhiVals(i),normStimCIEY(i)*max(priWhite_Y),'marker', '*','color',[nRGB0.red nRGB0.green nRGB0.blue]);    
    xlim([-5 350]); xlabel('phi'); ylabel('Y');ylim([0 120]);
    title('CIE Y across Azimuth (DKL)') ;   
            
    
    % Compare Lablib and Matlab colors. Directly plot the rgb coordinates
    subplot(235);
    rgb.red = rgb1(i,1); rgb.green = rgb1(i,2); rgb.blue = rgb1(i,3);
    nRGB1 = normalizeColors(rgb,kdlPhi,kdlTheta);
    
    plot3(rgb0(i,1),rgb0(i,2),rgb0(i,3),'marker', '+','color',[nRGB0.red nRGB0.green nRGB0.blue]);
    hold on;
    plot3(rgb1(i,1),rgb1(i,2),rgb1(i,3),'marker', 'o','color',[nRGB1.red nRGB1.green nRGB1.blue]);
    xlabel('r'); ylabel('g'); zlabel('b');
    xlim([0 1]); ylim([0 1]); zlim([0 1]);
    title('RGB values. Lablib (+) and Matlab (o)');
    
    subplot(236);
    % xyY coordinates projected back from the rgb space
    plot3(x,y,0.5+lum/2,'marker', '.','color',[nRGB0.red nRGB0.green nRGB0.blue]);
    hold on;

    % from Lablib RGB
    XYZlablib = Md65*[nRGB0.red nRGB0.green nRGB0.blue]'; xyYlablib = XYZToxyY(XYZlablib);
    plot3(xyYlablib(1),xyYlablib(2),xyYlablib(3),'marker', '+','color',[nRGB0.red nRGB0.green nRGB0.blue]);
    
    % from Matlab RGB
    XYZMatlab = Md65*[nRGB1.red nRGB1.green nRGB1.blue]'; xyYMatlab = XYZToxyY(XYZMatlab);
    plot3(xyYMatlab(1),xyYMatlab(2),xyYMatlab(3),'marker', 'o','color',[nRGB1.red nRGB1.green nRGB1.blue]);
    title('xy values from RGB. Desired (.), Lablib (+), Matlab (o)');
    xlabel('x'); ylabel('y'); zlabel('Y');
    xlim([0 1]); ylim([0 1]); zlim([0 1]);
end

end

% Macleod Boynton
function CMF2CF_MB = getMBMatrix
CMF2CF_MB = [0.15514 0.54312  -0.03286;
    -0.15514 0.45684  0.03286;
    0       0        0.01608];
end
function [r,b] = xy2MB(x,y)

CMF2CF_MB = getMBMatrix;
lms = CMF2CF_MB*[x y 1-x-y]';
r = lms(1)/(lms(1)+lms(2));
b = lms(3)/(lms(1)+lms(2));
end
function [x,y] = MB2xy(r,b)

CMF2CF_MB = getMBMatrix;
XYZ = (CMF2CF_MB)\[r 1-r b]';
sumXYZ = sum(XYZ);
x = XYZ(1)/sumXYZ;
y = XYZ(2)/sumXYZ;

end
% Lablib functions
function calibratedColor = computeKDLColors(bmr,bmb)

% cardinal green
cardg = (bmb(4) - bmb(3)) / (bmb(2) - bmb(3));

calibratedColor.cardinalGreen.red = 0.0;
calibratedColor.cardinalGreen.green = cardg;
calibratedColor.cardinalGreen.blue = 1 - cardg;

% cardinal yellow
cardy = (bmr(4) - bmr(2)) / (bmr(1) - bmr(2));

calibratedColor.cardinalYellow.red = cardy;
calibratedColor.cardinalYellow.green = 1 - cardy;
calibratedColor.cardinalYellow.blue = 0.0;

% equal energy

gbrb = (bmb(2) - bmb(3)) / (bmr(2) - bmr(3));

rede = bmb(4) - bmb(3) - (bmr(4) - bmr(3)) * gbrb;
redeDeno = bmb(1) - bmb(3) - (bmr(1) - bmr(3)) * gbrb;
rede = rede/redeDeno;

greene = (bmb(3) - bmb(4)) / (bmb(3) - bmb(2)) + rede * (bmb(1) - bmb(3)) / (bmb(3) - bmb(2));
bluee = 1 - rede - greene;

calibratedColor.equalEnergy.red = rede;
calibratedColor.equalEnergy.green = greene;
calibratedColor.equalEnergy.blue = bluee;

end
function kdlConstants = calibratedColorTokdlConstants(calibratedColor)

kdlConstants.rtc = (calibratedColor.equalEnergy.red - calibratedColor.cardinalYellow.red) /...
    calibratedColor.equalEnergy.red;
kdlConstants.gcb = (calibratedColor.equalEnergy.green - calibratedColor.cardinalGreen.green) /...
    calibratedColor.equalEnergy.green;
kdlConstants.gtc = (calibratedColor.equalEnergy.green - calibratedColor.cardinalYellow.green) /...
    calibratedColor.equalEnergy.green;
kdlConstants.bcb = (calibratedColor.equalEnergy.blue - calibratedColor.cardinalGreen.blue) /...
    calibratedColor.equalEnergy.blue;

end
function [rgb,lum,cb,tc] = kdlToRGB(kdlConstants,kdlPhi,kdlTheta)

% normFactor = sqrt(2.0);
normFactor = 1;
lum = sin(deg2rad(kdlTheta))/normFactor;
cb  = cos(deg2rad(kdlTheta)); tc=cb;
cb  = cb*cos(deg2rad(kdlPhi))/normFactor;
tc  = tc*sin(deg2rad(kdlPhi))/normFactor;

rgb.red   = (lum + cb + kdlConstants.rtc * tc);
rgb.green = (lum + cb * kdlConstants.gcb + tc * kdlConstants.gtc);
rgb.blue  = (lum + cb * kdlConstants.bcb + tc);

%normalize the values between 0 and 1;

% rgb.red = (rgb.red+1)/2;
% rgb.green = (rgb.green+1)/2;
% rgb.blue  = (rgb.blue+1)/2;

% original normalization followed inside the lablib for rendering using opengl

rgb.red = 0.5-(rgb.red/2);
rgb.green = 0.5-(rgb.green/2);
rgb.blue = 0.5-(rgb.blue/2);
end
function rgb = normalizeColors(rgb,kdlPhi,kdlTheta)

if (rgb.red>1) || (rgb.red<0)
    x = rgb.red;
    x = min(max(x,0),1);
    disp(['for phi=' num2str(kdlPhi) ', theta=' num2str(kdlTheta) ', red: ' num2str(rgb.red) ' out of range, set to ' num2str(x)]);
    rgb.red=x;
end
if (rgb.green>1) || (rgb.green<0)
    x = rgb.green;
    x = min(max(x,0),1);
    disp(['for phi=' num2str(kdlPhi) ', theta=' num2str(kdlTheta) ', green: ' num2str(rgb.green) ' out of range, set to ' num2str(x)]);
    rgb.green=x;
end
if (rgb.blue>1) || (rgb.blue<0)
    x = rgb.blue;
    x = min(max(x,0),1);
    disp(['for phi=' num2str(kdlPhi) ', theta=' num2str(kdlTheta) ', blue: ' num2str(rgb.blue) ' out of range, set to ' num2str(x)]);
    rgb.blue=x;
end

end
