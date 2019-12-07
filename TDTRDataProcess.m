function [] = TDTRDataProcess(varargin)
% Main function for the TDTR data processing.
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %
% TDTRDataProcess
% Main function for the TDTR data processing.
% This function can get three input, the dir for configuration file, the
% dir for data file/folder, the dir for output folder
% Reference:
% 1. main_TDTR_171120 in TDTR_Iwamoto_171120
% 2. The thesis of Aaron Jerome Schmidt
% Author: RL
% Date: Nov. 12, 2019
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %
close all
withInput = nargin;

% display some information about the code
disp('Thanks for using this code to process your data!');
disp('This code has several modes:')
disp('1. Fitting');
disp('2. Theory Curve');
disp('3. Sensitivity');
disp('4. Two-Frequency Fitting');
disp('5. Uncertainty');
disp('Note: If some errors existing, the file made by this code can not be deleted, you can run <fclose all> to solve this problem.');

global config
set(0,'defaultAxesFontSize',30)
set(0,'defaultTextFontSize',30)
set(0,'defaultAxesFontName','Helvetica')
set(0,'defaultTextFontName','Helvetica')
fPosition = [1 1 1400 1200];
% the configuration file
if withInput >= 1
    config_file = varargin{1};
    [path, filename, ext] = fileparts(config_file);
    filename = [filename, ext];
else
    [filename, path] = uigetfile('*.tdtrcfg','Pick a configuration file');
    if isequal(filename,0) || isequal(path,0)
        disp('User pressed cancel')
        return
    end
    config_file = fullfile(path, filename);
end
% a file with a suffix of .tdtrcfg can't be run directly, so we need to
% change its suffix to .m
temp_path = fullfile(path, 'temp');
if exist(temp_path,'dir') == 0
    mkdir(temp_path);
end
config_file_m = fullfile(temp_path, [filename(1:end-7), 'm']);
    copyfile(config_file, config_file_m);
config.fitting_mode = 0;
config.TheoryCurve_mode = 0;
config.Sensitivity = 0;
config.TwoFrequencyFitting = 0;
config.Uncertainty = 0;
run(config_file_m);

% fitting the data to get wanted parameters
if config.fitting_mode == 1
    length_para = size(config.fit_para,1);
    [para_name, format_f] = getLabel(config.LayerName, config.fit_para);
    % all data files in the selected folder will be processed
    if config.folder_mode == 1
        % the folder or file path of the data
        if withInput >= 2
            SourcePath = varargin{2};
        else
            SourcePath = uigetdir('.','Pick a folder where there are data files');
            if isequal(SourcePath,0)
                disp('User pressed cancel')
                return
            end
        end
        % the folder to storage the result files
        if withInput >= 3
            OutputPath = varargin{3};
        else
            OutputPath = uigetdir('.','Pick a folder to storage the results');
            if isequal(OutputPath,0)
                disp('User pressed cancel')
                return
            end
        end
        IsNotExist = 0;
        index = 1;
        while IsNotExist == 0
            OutputFolder = [datestr(datetime('now'),'yyyy-mm-dd_HH-MM-SS'), '__', num2str(index)];
            OutputFolder = fullfile(OutputPath, OutputFolder);
            if exist(OutputFolder,'dir') == 0
                IsNotExist = 1;
            else
                index = index + 1;
            end
        end
        % get the files in the SourcePath
        filelist = dir([SourcePath '\*.txt']);
        length_filelist = length(filelist);
        if length_filelist < 1
            disp('There is not data file in selected folder.');
            return
        end
        Results = zeros(length_filelist,length_para);
        % make output folder
        % make the root folder for output
        mkdir(OutputFolder);
        % make the image folder
        img_folder = fullfile(OutputFolder,'img');
        mkdir(img_folder);
        % make the raw data folder
        % there are three columns in each file, which are time(ns), x, y
        raw_data_folder = fullfile(OutputFolder,'raw_data');
        mkdir(raw_data_folder);
        % make the dealed raw data folder
        % there are two columns on each file, which are time(ns), fun
        dealed_data_folder = fullfile(OutputFolder,'dealed_data');
        mkdir(dealed_data_folder);
        % the theory data folder
        % there are two columns on each file, which are time(ns), fun
        theory_data_folder = fullfile(OutputFolder,'theory_data');
        mkdir(theory_data_folder);
        % copy the configuration file to output folder
        copyfile(config_file, fullfile(OutputFolder,'configuration.txt'));
        % make result_log file in output folder to save fitting results
        result_file = fopen(fullfile(OutputFolder,'result_log.txt'),'a+');
        
        disp(['There are ', num2str(length_filelist), ' data files will be processed.']);
        fprintf(result_file,'%s\r\n',['There are ', num2str(length_filelist), ' data files will be processed.']);
        for index = 1:1:length_filelist
            disp(['File ', num2str(index), ' : ', filelist(index).name]);
            fprintf(result_file,'%5s %3s %3s %s\r\n','File', num2str(index), ':', filelist(index).name);
            % read raw data
            raw_data = load(fullfile(SourcePath, filelist(index).name));
            % do the fitting and get result
            Result = TDTRDataFitting(raw_data, config);
            % save the raw data file to output folder
            raw_data_file = fullfile(raw_data_folder,filelist(index).name);
            copyfile(fullfile(SourcePath, filelist(index).name), raw_data_file);
            % save the dealed data file to output folder
            dealed_data_file = fopen(fullfile(dealed_data_folder, filelist(index).name),'w+');
            fprintf(dealed_data_file, '%f\t%f\r\n', [Result.dealed_data.tau(:)'; Result.dealed_data.fun(:)']);
            fclose(dealed_data_file);
            % save the theory data to output file
            theory_data_file = fopen(fullfile(theory_data_folder, filelist(index).name),'w+');
            fprintf(theory_data_file, '%f\t%f\r\n', [Result.theory_data.tau(:)'; Result.theory_data.fun(:)']);
            fclose(theory_data_file);
            % save the figure of dealed raw data and theory data to output file
            fig = figure('Position', fPosition);
            semilogx(Result.dealed_data.tau,Result.dealed_data.fun,'ko','MarkerSize',8);
            title(['Fitting: ', num2str(index)]);
            hold on
            semilogx(Result.theory_data.tau,Result.theory_data.fun,'r-','LineWidth',4);
            if config.mode ~= 'p'
                ylim([0,max(Result.dealed_data.fun)+0.5]) 
            end
            legend('Data','Best fit')
            xlabel('Delay time [s]')
            switch config.mode
                case 'r'
                    ylabel('Ratio')
                case 'p'
                    ylabel('Phase')
                case 'a'
                    ylabel('Amplitude')
            end
            saveas(fig,fullfile(img_folder,[filelist(index).name(1:end-4),'.png']),'png');
            % diaplay the result of fitting
            disp(para_name);
            disp(Result.fittingValue);
            disp(['The standard deviation is ', num2str(Result.StdDev)]);
            % save the fitting result to result_log file
            fprintf(result_file,'%s\r\n',para_name);
            fprintf(result_file,format_f,Result.fittingValue');
            fprintf(result_file,'%s\r\n',['The standard deviation is ', num2str(Result.StdDev)]);
            Results(index,:) = Result.fittingValue;
        end
        disp('Summary of the results');
        disp(para_name);
        disp(Results);
        disp('The average value');
        disp(mean(Results,1));
        disp('The standard deviation');
        disp(std(Results,0,1));
        fprintf(result_file,'\r\n%s\r\n','Summary of the results');
        fprintf(result_file,'%s\r\n',para_name);
        fprintf(result_file,format_f,Results');
        fprintf(result_file,'%s\r\n','The average value');
        fprintf(result_file,format_f,mean(Results,1));
        fprintf(result_file,'%s\r\n','The standard deviation');
        fprintf(result_file,format_f,std(Results,0,1));
        fclose(result_file);
        
    % a selected file will be processed
    else
        % the folder or file path of the data
        if withInput >= 2
            SourcePath = varargin{2};
        else
            [filename, path] = uigetfile('*.txt','Pick a data file');
            if isequal(filename,0) || isequal(path,0)
                disp('User pressed cancel')
                return
            end
            SourcePath = fullfile(path, filename);
        end
        % the folder to storage the result files
        if withInput >= 3
            OutputPath = varargin{3};
        else
            OutputPath = uigetdir('.','Pick a folder to storage the results');
            if isequal(OutputPath,0)
                disp('User pressed cancel')
                return
            end
        end
        IsNotExist = 0;
        index = 1;
        while IsNotExist == 0
            OutputFolder = [datestr(datetime('now'),'yyyy-mm-dd_HH-MM-SS'), '__', num2str(index)];
            OutputFolder = fullfile(OutputPath, OutputFolder);
            if exist(OutputFolder,'dir') == 0
                IsNotExist = 1;
            else
                index = index + 1;
            end
        end        
        % make output folder
        % Different from folder mode, the result files will be put in output folder directly in single file mode 
        mkdir(OutputFolder);
        % copy the configuration file to output folder
        copyfile(config_file, fullfile(OutputFolder,'configuration.txt'));
        % make result_log file in output folder to save fitting results
        result_file = fopen(fullfile(OutputFolder,'result_log.txt'),'a+');
                
        disp(['Filename: ', filename]);
        fprintf(result_file,'Filename: %s\r\n',filename);
        % read raw data
        raw_data = load(SourcePath);
        % do the fitting and get result
        Result = TDTRDataFitting(raw_data, config);
        % save the raw data file to output folder
        raw_data_file = fullfile(OutputFolder,[filename(1:end-4), '_raw_data.txt']);
        copyfile(SourcePath, raw_data_file);
        % save the dealed data file to output folder
        dealed_data_file = fopen(fullfile(OutputFolder, [filename(1:end-4), '_dealed_data.txt']),'w+');
        fprintf(dealed_data_file, '%f\t%f\r\n', [Result.dealed_data.tau(:)'; Result.dealed_data.fun(:)']);
        fclose(dealed_data_file);
        % save the theory data to output file
        theory_data_file = fopen(fullfile(OutputFolder, [filename(1:end-4), '_theory_data.txt']),'w+');
        fprintf(theory_data_file, '%f\t%f\r\n', [Result.theory_data.tau(:)'; Result.theory_data.fun(:)']);
        fclose(theory_data_file);
        % save the figure of dealed raw data and theory data to output file
        fig = figure('Position', fPosition);
        semilogx(Result.dealed_data.tau,Result.dealed_data.fun,'ko','MarkerSize',8);
        title(['Fitting: ', num2str(index)]);
        hold on
        semilogx(Result.theory_data.tau,Result.theory_data.fun,'r-','LineWidth',4);
        if config.mode ~= 'p'
            ylim([0,max(Result.dealed_data.fun)+0.5]) 
        end
        legend('Data','Best fit')
        xlabel('Delay time [s]')
        switch config.mode
            case 'r'
                ylabel('Ratio')
            case 'p'
                ylabel('Phase')
            case 'a'
                ylabel('Amplitude')
        end
        saveas(fig,fullfile(OutputFolder,[filename(1:end-4),'.png']),'png');
        % diaplay the result of fitting
        disp(para_name);
        disp(Result.fittingValue);
        disp(['The standard deviation is ', num2str(Result.StdDev)]);
        % save the fitting result to result_log file
        fprintf(result_file,'%s\r\n',para_name);
        fprintf(result_file,format_f,Result.fittingValue);
        fprintf(result_file,'%s\r\n',['The standard deviation is ', num2str(Result.StdDev)]);
        fclose(result_file);
    end   
end
%
% Plot the theory curve
if config.TheoryCurve_mode == 1
    tau_data = logspace(log10(config.tau(1)),log10(config.tau(2)),200);
    [func] = TheoryData(config.kz,config.kr,config.G,config.d,config.vhc,config.w,tau_data, config);
    fig = figure('Position', fPosition);
    semilogx(tau_data,func,'ko','MarkerSize',8);
    title('Theory curve');
    OutputPath = uigetdir('.','Pick a folder to storage the results');
    if isequal(OutputPath,0)
        disp('User pressed cancel')
        return
    end
    IsNotExist = 0;
    index = 1;
    while IsNotExist == 0
        OutputFolder = [datestr(datetime('now'),'yyyy-mm-dd_HH-MM-SS'), '__', num2str(index)];
        OutputFolder = fullfile(OutputPath, OutputFolder);
        if exist(OutputFolder,'dir') == 0
            IsNotExist = 1;
        else
            index = index + 1;
        end
    end
    mkdir(OutputFolder);
    copyfile(config_file, fullfile(OutputFolder,'configuration.txt'));
    % save the theory data to output file
    theory_data_file = fopen(fullfile(OutputFolder, 'theory_data.txt'),'w+');
    fprintf(theory_data_file, '%f\t%f\r\n', [tau_data(:)'*1E9; func(:)']);
    fclose(theory_data_file);
end

if config.TwoFrequencyFitting == 1
    Ndata = length(config.Data);
    length_para = zeros(1, Ndata);
    for index = 1:1:Ndata
        [para_name{index}, format_f{index}] = getLabel(config.Data{index}.LayerName, config.Data{index}.fit_para);
        length_para(index) = size(config.Data{index}.fit_para,1);
    end
    SourcePath_set = cell(1,Ndata);
    % all data files in the selected folder will be processed
    if config.folder_mode == 1
        if (withInput > 1)
            if (withInput < (Ndata+2))
                disp('The number of input parameters is not enough');
                return
            else
                for index = 1:1:Ndata
                    SourcePath_set{index} = varargin{index+1};
                end
                OutputPath = varargin{index+2};
            end
        else
            for index = 1:1:Ndata
                SourcePath_set{index} = uigetdir('.',['Pick a folder where there are data files for ', num2str(index),'th data.']);
                if isequal(SourcePath_set{index},0)
                    disp('User pressed cancel')
                    return
                end
            end
            OutputPath = uigetdir('.','Pick a folder to storage the results');
            if isequal(OutputPath,0)
                disp('User pressed cancel')
                return
            end
        end
        IsNotExist = 0;
        index = 1;
        while IsNotExist == 0
            OutputFolder = [datestr(datetime('now'),'yyyy-mm-dd_HH-MM-SS'), '__', num2str(index)];
            OutputFolder = fullfile(OutputPath, OutputFolder);
            if exist(OutputFolder,'dir') == 0
                IsNotExist = 1;
            else
                index = index + 1;
            end
        end
        filelist_set = cell(1,Ndata);
        length_filelist_set = zeros(1, Ndata);
        for index = 1:1:Ndata
            filelist_set{index}=dir(fullfile(SourcePath_set{index}, '*.txt'));
            length_filelist_set(index)=length(filelist_set{index});
            if length_filelist_set(index) < 1
                disp('There is not data file in selected folder.');
                return
            end
        end
        length_filelist = min(length_filelist_set);
        Result_set = cell(1,Ndata);
        for index = 1:1:Ndata
            Result_set{index} = zeros(length_filelist,length_para(index));
        end
        % make output folder
        % make the root folder for output
        mkdir(OutputFolder);
        OutputFolder_set = cell(1,Ndata);
        img_folder_set = cell(1,Ndata);
        raw_data_folder_set = cell(1,Ndata);
        dealed_data_folder_set = cell(1,Ndata);
        theory_data_folder_set = cell(1,Ndata);
        result_file_set = cell(1,Ndata);
        for index = 1:1:Ndata
            OutputFolder_set{index} = fullfile(OutputFolder, ['Data_1', num2str(index)]);
            % make the image folder
            img_folder_set{index} = fullfile(OutputFolder_set{index},'img');
            mkdir(img_folder_set{index});
            % make the raw data folder
            % there are three columns in each file, which are time(ns), x, y
            raw_data_folder_set{index} = fullfile(OutputFolder_set{index},'raw_data');
            mkdir(raw_data_folder_set{index});
            % make the dealed raw data folder
            % there are two columns on each file, which are time(ns), fun
            dealed_data_folder_set{index} = fullfile(OutputFolder_set{index},'dealed_data');
            mkdir(dealed_data_folder_set{index});
            % the theory data folder
            % there are two columns on each file, which are time(ns), fun
            theory_data_folder_set{index} = fullfile(OutputFolder_set{index},'theory_data');
            mkdir(theory_data_folder_set{index});
            % make result_log file in output folder to save fitting results
            result_file_set{index} = fopen(fullfile(OutputFolder_set{index},'result_log.txt'),'a+');
        end   
        % copy the configuration file to output folder
        copyfile(config_file, fullfile(OutputFolder,'configuration.txt'));
        disp(['There are ', num2str(length_filelist), ' data files will be processed.']);
        for index = 1:1:Ndata
            fprintf(result_file_set{index},'%s\r\n',['There are ', num2str(length_filelist), ' data files will be processed.']);
        end
        for index = 1:1:length_filelist
            disp(['File ', num2str(index), ' : ']);
            raw_data_set = cell(1, Ndata);
            for index_1 = 1:1:Ndata
                fprintf(result_file_set{index_1},'%5s %3s %3s %s\r\n','File', num2str(index_1), ':', filelist_set{index_1}(index).name);
                % read raw data
                raw_data_set{index_1} = load(fullfile(SourcePath_set{index_1}, filelist_set{index_1}(index).name));
            end
            % do the fitting and get result
            Result = TDTRDataTwoFreFitting(raw_data_set, config);
            for index_1 = 1:1:Ndata
                % save the raw data file to output folder
                raw_data_file = fullfile(raw_data_folder_set{index_1},filelist_set{index_1}(index).name);
                copyfile(fullfile(SourcePath_set{index_1}, filelist_set{index_1}(index).name), raw_data_file);
                % save the dealed data file to output folder
                dealed_data_file = fopen(fullfile(dealed_data_folder_set{index_1}, filelist_set{index_1}(index).name),'w+');
                fprintf(dealed_data_file, '%f\t%f\r\n', [Result{index_1}.dealed_data.tau(:)'; Result{index_1}.dealed_data.fun(:)']);
                fclose(dealed_data_file);
                % save the theory data to output file
                theory_data_file = fopen(fullfile(theory_data_folder_set{index_1}, filelist_set{index_1}(index).name),'w+');
                fprintf(theory_data_file, '%f\t%f\r\n', [Result{index_1}.theory_data.tau(:)'; Result{index_1}.theory_data.fun(:)']);
                fclose(theory_data_file);
                % save the figure of dealed raw data and theory data to output file
                fig = figure('Position', fPosition);
                semilogx(Result{index_1}.dealed_data.tau,Result{index_1}.dealed_data.fun,'ko','MarkerSize',8);
                title(['Fitting: ', num2str(index)]);
                hold on
                semilogx(Result{index_1}.theory_data.tau,Result{index_1}.theory_data.fun,'r-','LineWidth',4);
                if config.Data{index_1}.mode ~= 'p'
                    ylim([0,max(Result{index_1}.dealed_data.fun)+0.5]) 
                end
                legend('Data','Best fit')
                xlabel('Delay time [s]')
                switch config.Data{index_1}.mode
                    case 'r'
                        ylabel('Ratio')
                    case 'p'
                        ylabel('Phase')
                    case 'a'
                        ylabel('Amplitude')
                end
                saveas(fig,fullfile(img_folder_set{index_1},[filelist_set{index_1}(index).name(1:end-4),'.png']),'png');
                % diaplay the result of fitting
                disp(para_name{index_1});
                disp(Result{index_1}.fittingValue);
                disp(['The standard deviation is ', num2str(Result{index_1}.StdDev)]);
                % save the fitting result to result_log file
                fprintf(result_file_set{index_1},'%s\r\n',para_name{index_1});
                fprintf(result_file_set{index_1},format_f{index_1},Result{index_1}.fittingValue');
                fprintf(result_file_set{index_1},'%s\r\n',['The standard deviation is ', num2str(Result{index_1}.StdDev)]);
                Result_set{index_1}(index,:) = Result{index_1}.fittingValue;
            end
        end
        disp('Summary of the results');
        for index = 1:1:Ndata
            disp(['Data',num2str(index),':']);
            disp(para_name{index});
            disp(Result_set{index});
            disp('The average value');
            disp(mean(Result_set{index},1));
            disp('The standard deviation');
            disp(std(Result_set{index},0,1));
            fprintf(result_file_set{index},'\r\n%s\r\n','Summary of the results');
            fprintf(result_file_set{index},'%s\r\n',para_name{index});
            fprintf(result_file_set{index},format_f{index},Result_set{index}');
            fprintf(result_file_set{index},'%s\r\n','The average value');
            fprintf(result_file_set{index},format_f{index},mean(Result_set{index},1));
            fprintf(result_file_set{index},'%s\r\n','The standard deviation');
            fprintf(result_file_set{index},format_f{index},std(Result_set{index},0,1));
            fclose(result_file_set{index});
        end
    else
        
    end
end

% close file and delete the temp folder
fclose all;
rmdir(temp_path, 's');