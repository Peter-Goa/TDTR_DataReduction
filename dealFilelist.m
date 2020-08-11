function [X_max, Y_max, Index_list] = dealFilelist(filelist)

length_filelist = length(filelist);
Index_list = zeros(length_filelist,2);
for index = 1:1:length_filelist
    filename = filelist(index).name;
    exp1 = '\(\d+,\d+\)';
    exp2 = '\d+';
    [start_index, end_index] = regexp(filename,exp1);
    temp_str = filename(start_index:end_index);
    [start_index, end_index] = regexp(temp_str,exp2);
    Index_list(index,1) = str2double(temp_str(start_index(1):end_index(1)))+1;
    Index_list(index,2) = str2double(temp_str(start_index(2):end_index(2)))+1;
end
X_max = max(Index_list(:,1));
Y_max = max(Index_list(:,2));

