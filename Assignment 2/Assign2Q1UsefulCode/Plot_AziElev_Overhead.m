% Plot_AziElev_Overhead - Script file to setup a MATLAB figure for simple
% plotting of an azimuth-elevation overhead plot. Also use
% "overheadplotcoords" to transform azimuth and elevation angles into plot
% coordinates for this figure 

% Setup Figure
figure
hold on

% Plot Markings on Figure
tincs = 0:0.2:(2*pi + 0.2);
xpol(1,:) = 30*sin(tincs);
xpol(2,:) = 60*sin(tincs);
xpol(3,:) = 90*sin(tincs);
ypol(1,:) = 30*cos(tincs);
ypol(2,:) = 60*cos(tincs);
ypol(3,:) = 90*cos(tincs);
plot([0,0],[-90,90])
plot([-90,90],[0,0])
plot(xpol(1,:),ypol(1,:))
plot(xpol(2,:),ypol(2,:))
plot(xpol(3,:),ypol(3,:))

% Plot Directions and Numbers on Figures
text(0,90,'North')
text(-90,0,'East')
text(0,-90,'South')
text(90,0,'West')
text(0,0,'90^o')
text(30*0.7,30*0.7,'60^o')
text(60*0.7,60*0.7,'30^o')
text(90*0.7,90*0.7,'0^o')

