% This code is used to plot the data of TDTR measurement. It can output
% A-t (not normalized)and R-t picture.
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
end
if exist([filepath,'\img_R'],'dir')==0
    mkdir(filepath,'img_R');
end
Color_list = ['b','r','m','k'];
Dot_list = ['o','x','+','*','s','d','v','^','<','>','p','h'];
legend_list = cell(length_namelist,1);
for index = 1:length_namelist
    legend_list{index} = namelist(index).name(1:end-4);
end
figure(1)
clf(1)
hold on
for index = 1:length_namelist
    %disp(namelist(index).name);
    data = load([filepath '\' namelist(index).name]);
    %disp(data)
    Marker_type = [Dot_list(fix(index/4)+1), Color_list(mod(index,4)+1)];
    plot(data(:,1),sqrt(data(:,2).^2+data(:,3).^2),Marker_type);
end
legend(legend_list,'Location','BestOutside');
xlabel('Time(ns)');
ylabel('Amplitude');
set(gcf,'position',[300, 300,1000, 700]);
hold off
saveas(gcf,[filepath,'\img_A\', 'A_all', '.jpg']);

figure(2)
clf(2)
hold on
for index = 1:length_namelist
    %disp(namelist(index).name);
    data = load([filepath '\' namelist(index).name]);
    %disp(data)
    Marker_type = [Dot_list(fix(index/4)+1), Color_list(mod(index,4)+1)];
    plot(data(:,1),abs(data(:,2)./data(:,3)),Marker_type);
end
legend(legend_list,'Location','BestOutside');
xlabel('Time(ns)');
ylabel('Ratio');
set(gcf,'position',[300, 300,1000, 700]);
hold off
saveas(gcf,[filepath,'\img_R\', 'R_all', '.jpg']);
