% This code is used to plot the data of TDTR measurement. It can output
% all dara in one A-t (not normalized)and one R-t picture.
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %
% Plot TDTR Data
% Author: Rulei
% Date: Oct. 23th 2019
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %

% select a direction by UI
filepath = uigetdir('*.*','select a direction');
namelist = dir([filepath '\*.txt']);
length_namelist = length(namelist);
if exist([filepath,'\img_A'],'dir')==0
    mkdir(filepath,'img_A');
else
    disp('Dirction img_A exists, so can not create this direction again.');
%   return
end
if exist([filepath,'\img_R'],'dir')==0
    mkdir(filepath,'img_R');
else
    disp('Dirction img_R exists, so can not create this direction again.');
%   return
end

for index = 1:length_namelist
    %disp(namelist(index).name);
    data = load([filepath '\' namelist(index).name]);
    %disp(data)
    plot(data(:,1),sqrt(data(:,2).^2+data(:,3).^2),'oR');
    xlabel('Time(ns)');
    ylabel('Amplitude');
    set(gcf,'position',[300, 300,1000, 700])
    saveas(gcf,[filepath,'\img_A\', namelist(index).name(1:end-4), '.jpg']);
    plot(data(:,1),abs(data(:,2)./data(:,3)),'*b');
    xlabel('Time(ns)');
    ylabel('Ratio');
    saveas(gcf,[filepath,'\img_R\', namelist(index).name(1:end-4), '.jpg']);
end