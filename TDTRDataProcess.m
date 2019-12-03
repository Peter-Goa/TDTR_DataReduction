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
else
    [filename, path] = uigetfile('*.tdtrcfg','Pick a configuration file');
    if isequal(filename,0) || isequal(path,0)
        disp('User pressed cancel')
        return
    end
    config_file = fullfile(path, filename);
    % a file with a suffix of .tdtrcfg can't be run directly, so we need to
    % change its suffix to .m
    temp_path = fullfile(path, 'temp');
    if exist(temp_path,'dir') == 0
        mkdir(temp_path);
    end
    config_file_m = fullfile(temp_path, [filename(1:end-7), 'm']);
    copyfile(config_file, config_file_m);
end
config.fitting_mode = 0;
config.TheoryCurve_mode = 0;
config.Sensitivity = 0;
config.TwoFrequencyFitting = 0;
config.Uncertainty = 0;
run(config_file_m);

% fitting the data to get wanted parameters
if config.fitting_mode == 1
    length_para = size(config.fit_para,1);
    para_name = '';
    for index_para = 1:1:length_para
        para_name = [para_name, config.LayerName{config.fit_para(index_para,1)}, '_'];
        switch config.fit_para(index_para, 2)
            case 1
                para_name = [para_name, 'kz [W/mK]'];
            case 2
                para_name = [para_name, 'kr [W/mK]'];
            case 3
                para_name = [para_name, 'vhc [MJ/m^3K]'];
            case 4
                para_name = [para_name, 'd [nm]'];
            case 5
                para_name = [para_name, 'G [MW/m^2K]'];
        end
        para_name = [para_name, '   ']; 
    end

    format_f = '';
    for index = 1:1:length_para
        format_f = [format_f '%f  '];
    end
    format_f = [format_f '\r\n'];
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
            Result = TDTRDataFitting(raw_data);
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
            saveas(fig,fullfile(img_folder,filelist(index).name(1:end-4)),'png');
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
        Result = TDTRDataFitting(raw_data);
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
        saveas(fig,fullfile(OutputFolder,filename(1:end-4)),'png');
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

% close file and delete the temp folder
fclose all;
rmdir(temp_path, 's');