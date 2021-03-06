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
disp('6. Electron-Phonon Coupling factor');
disp('7. Mapping');
disp('Note: If some errors existing, the file made by this code can not be deleted, you can run <fclose all> to solve this problem.');

global config;
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
config.ZeroPointMode = 0;
config.ElectronPhononMode = 0;
config.mapping_mode = 0;
config.ParameterScan = 0;
run(config_file_m);

%% fitting the data to get wanted parameters
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
                fclose all;
                rmdir(temp_path, 's');
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
                fclose all;
                rmdir(temp_path, 's');
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
        filelist = dir([SourcePath filesep() '*.txt']);
        length_filelist = length(filelist);
        if length_filelist < 1
            disp('There is not data file in selected folder.');
            fclose all;
            rmdir(temp_path, 's');
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
                fclose all;
                rmdir(temp_path, 's');
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
                fclose all;
                rmdir(temp_path, 's');
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

%% Plot the theory curve
if config.TheoryCurve_mode == 1
    tau_data = logspace(log10(config.tau(1)),log10(config.tau(2)),200);
    [func] = TheoryData(config.kz,config.kr,config.G,config.d,config.vhc,config.w,tau_data, config);
    fig = figure('Position', fPosition);
    semilogx(tau_data,func,'ko','MarkerSize',8);
    title('Theory curve');
    if config.Bias == 1
        hold on;
        switch config.Bias_para(2)
            case 1
                kz_temp_plus = config.kz;
                kz_temp_plus(config.Bias_para(1)) = kz_temp_plus(config.Bias_para(1))*(1+config.Bias_para(3)/100);
                [func_plus] = TheoryData(kz_temp_plus,config.kr,config.G,config.d,config.vhc,config.w,tau_data, config);
                kz_temp_minus = config.kz;
                kz_temp_minus(config.Bias_para(1)) = kz_temp_minus(config.Bias_para(1))*(1-config.Bias_para(3)/100);
                [func_minus] = TheoryData(kz_temp_minus,config.kr,config.G,config.d,config.vhc,config.w,tau_data, config);
                Label = {'origin',[num2str(100+config.Bias_para(3)),'% ', num2str(config.Bias_para(1)),'th kz'], [num2str(100-config.Bias_para(3)),'% ', num2str(config.Bias_para(1)),'th kz']};
            case 2
                kr_temp_plus = config.kr;
                kr_temp_plus(config.Bias_para(1)) = kr_temp_plus(config.Bias_para(1))*(1+config.Bias_para(3)/100);
                [func_plus] = TheoryData(config.kz,kr_temp_plus,config.G,config.d,config.vhc,config.w,tau_data, config);
                kr_temp_minus = config.kr;
                kr_temp_minus(config.Bias_para(1)) = kr_temp_minus(config.Bias_para(1))*(1-config.Bias_para(3)/100);
                [func_minus] = TheoryData(config.kz,kr_temp_minus,config.G,config.d,config.vhc,config.w,tau_data, config);
                Label = {'origin',[num2str(100+config.Bias_para(3)),'% ', num2str(config.Bias_para(1)),'th kr'], [num2str(100-config.Bias_para(3)),'% ', num2str(config.Bias_para(1)),'th kr']};
            case 3
                vhc_temp_plus = config.vhc;
                vhc_temp_plus(config.Bias_para(1)) = vhc_temp_plus(config.Bias_para(1))*(1+config.Bias_para(3)/100);
                [func_plus] = TheoryData(config.kz,config.kr,config.G,config.d,vhc_temp_plus,config.w,tau_data, config);
                vhc_temp_minus = config.vhc;
                vhc_temp_minus(config.Bias_para(1)) = vhc_temp_minus(config.Bias_para(1))*(1-config.Bias_para(3)/100);
                [func_minus] = TheoryData(config.kz,config.kr,config.G,config.d,vhc_temp_minus,config.w,tau_data, config);
                Label = {'origin',[num2str(100+config.Bias_para(3)),'% ', num2str(config.Bias_para(1)),'th vhc'], [num2str(100-config.Bias_para(3)),'% ', num2str(config.Bias_para(1)),'th vhc']};
            case 4
                d_temp_plus = config.d;
                d_temp_plus(config.Bias_para(1)) = d_temp_plus(config.Bias_para(1))*(1+config.Bias_para(3)/100);
                [func_plus] = TheoryData(config.kz,config.kr,config.G,d_temp_plus,config.vhc,config.w,tau_data, config);
                d_temp_minus = config.d;
                d_temp_minus(config.Bias_para(1)) = d_temp_minus(config.Bias_para(1))*(1-config.Bias_para(3)/100);
                [func_minus] = TheoryData(config.kz,config.kr,config.G,d_temp_minus,config.vhc,config.w,tau_data, config);
                Label = {'origin',[num2str(100+config.Bias_para(3)),'% ', num2str(config.Bias_para(1)),'th d'], [num2str(100-config.Bias_para(3)),'% ', num2str(config.Bias_para(1)),'th d']};                
            case 5
                G_temp_plus = config.G;
                G_temp_plus(config.Bias_para(1)) = G_temp_plus(config.Bias_para(1))*(1+config.Bias_para(3)/100);
                [func_plus] = TheoryData(config.kz,config.kr,G_temp_plus,config.d,config.vhc,config.w,tau_data, config);
                G_temp_minus = config.G;
                G_temp_minus(config.Bias_para(1)) = G_temp_minus(config.Bias_para(1))*(1-config.Bias_para(3)/100);
                [func_minus] = TheoryData(config.kz,config.kr,G_temp_minus,config.d,config.vhc,config.w,tau_data, config);
                Label = {'origin',[num2str(100+config.Bias_para(3)),'% ', num2str(config.Bias_para(1)),'th G'], [num2str(100-config.Bias_para(3)),'% ', num2str(config.Bias_para(1)),'th G']};       
        end
        semilogx(tau_data,func_plus,'r-.','LineWidth',2);
        semilogx(tau_data,func_minus,'b-.','LineWidth',2);
        legend(Label);
    end
    if withInput >= 2
        OutputPath = varargin{2};
    else
        OutputPath = uigetdir('.','Pick a folder to storage the results');
        if isequal(OutputPath,0)
            disp('User pressed cancel')
            fclose all;
            rmdir(temp_path, 's');
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
    mkdir(OutputFolder);
    copyfile(config_file, fullfile(OutputFolder,'configuration.txt'));
    % save the theory data to output file
    theory_data_file = fopen(fullfile(OutputFolder, 'theory_data.txt'),'w+');
    fprintf(theory_data_file, '%f\t%f\r\n', [tau_data(:)'*1E9; func(:)']);
    fclose(theory_data_file);
    hold off
    saveas(fig,fullfile(OutputFolder, 'TheoryCurve.png'),'png');
end

%% Two frequency fitting model
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
                fclose all;
                rmdir(temp_path, 's');
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
                    fclose all;
                    rmdir(temp_path, 's');
                    return
                end
            end
            OutputPath = uigetdir('.','Pick a folder to storage the results');
            if isequal(OutputPath,0)
                disp('User pressed cancel')
                fclose all;
                rmdir(temp_path, 's');
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
                fclose all;
                rmdir(temp_path, 's');
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
        filelist_set = cell(1,Ndata);
        if (withInput > 1)
            if (withInput < (Ndata+2))
                disp('The number of input parameters is not enough');
                fclose all;
                rmdir(temp_path, 's');
                return
            else
                for index = 1:1:Ndata
                    SourcePath_set{index} = varargin{index+1};
                end
                OutputPath = varargin{index+2};
            end
        else
            for index = 1:1:Ndata
                [filelist_set{index},SourcePath_set{index}] = uigetfile('*.txt',['Pick a data file for ', num2str(index),'th data.']);
                if isequal(SourcePath_set{index},0) || isequal(filelist_set{index},0)
                    disp('User pressed cancel')
                    fclose all;
                    rmdir(temp_path, 's');
                    return
                end
            end
            OutputPath = uigetdir('.','Pick a folder to storage the results');
            if isequal(OutputPath,0)
                disp('User pressed cancel')
                fclose all;
                rmdir(temp_path, 's');
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
        length_filelist = 1;
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
        raw_data_set = cell(1, Ndata);
        for index_1 = 1:1:Ndata
            fprintf(result_file_set{index_1},'%5s %3s %3s %s\r\n','File', num2str(index_1), ':', filelist_set{index_1});
            % read raw data
            raw_data_set{index_1} = load(fullfile(SourcePath_set{index_1}, filelist_set{index_1}));
        end
        % do the fitting and get result
        Result = TDTRDataTwoFreFitting(raw_data_set, config);
        for index_1 = 1:1:Ndata
            % save the raw data file to output folder
            raw_data_file = fullfile(raw_data_folder_set{index_1},filelist_set{index_1});
            copyfile(fullfile(SourcePath_set{index_1}, filelist_set{index_1}), raw_data_file);
            % save the dealed data file to output folder
            dealed_data_file = fopen(fullfile(dealed_data_folder_set{index_1}, filelist_set{index_1}),'w+');
            fprintf(dealed_data_file, '%f\t%f\r\n', [Result{index_1}.dealed_data.tau(:)'; Result{index_1}.dealed_data.fun(:)']);
            fclose(dealed_data_file);
            % save the theory data to output file
            theory_data_file = fopen(fullfile(theory_data_folder_set{index_1}, filelist_set{index_1}),'w+');
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
            saveas(fig,fullfile(img_folder_set{index_1},[filelist_set{index_1}(1:end-4),'.png']),'png');
            % diaplay the result of fitting
            disp(para_name{index_1});
            disp(Result{index_1}.fittingValue);
            disp(['The standard deviation is ', num2str(Result{index_1}.StdDev)]);
            % save the fitting result to result_log file
            fprintf(result_file_set{index_1},'%s\r\n',para_name{index_1});
            fprintf(result_file_set{index_1},format_f{index_1},Result{index_1}.fittingValue');
            fprintf(result_file_set{index_1},'%s\r\n',['The standard deviation is ', num2str(Result{index_1}.StdDev)]);
            Result_set{index_1} = Result{index_1}.fittingValue;
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
    end
end

%% Sensitivity calculation model
if config.Sensitivity == 1
    NumParame = size(config.sensitivity_para, 1);
    Numtime = 200;
    time = linspace(config.tau(1), config.tau(2), Numtime);
    time_ns = time*1e9;
    result = zeros(NumParame, Numtime);
    OriginalFun = TheoryData(config.kz,config.kr,config.G,config.d,config.vhc,config.w, time, config);
    para_name = getLabel_s(config.LayerName, config.sensitivity_para);
    para_name_s = 'Time[ns]';
    form_s = '%f\t';
    for index = 1:1:NumParame
        para_name_s = [para_name_s,',',para_name{index}];
        form_s = [form_s,'%f\t'];
    end
    form_s = [form_s,'\r\n'];
    Colors = getColors(NumParame);
    ratio = 1.01;
    Sens_Data = zeros(Numtime,NumParame+1);
    Sens_Data(:,1) = time_ns';
    for index = 1:1:NumParame
        kz = config.kz;
        kr = config.kr;
        G = config.G;
        d = config.d;
        vhc = config.vhc;
        w = config.w;
        if config.sensitivity_para(index,1) == 0
            switch config.sensitivity_para(index, 2)
                case 1
                    w(1) = w(1)*ratio;
                case 2
                    w(2) = w(2)*ratio;
            end
        else
            switch config.sensitivity_para(index, 2)
                case 1
                    kz(config.sensitivity_para(index,1)) = kz(config.sensitivity_para(index,1))*ratio;
                case 2
                    kr(config.sensitivity_para(index,1)) = kr(config.sensitivity_para(index,1))*ratio;
                case 3
                    vhc(config.sensitivity_para(index,1)) = vhc(config.sensitivity_para(index,1))*ratio;
                case 4
                    d(config.sensitivity_para(index,1)) = d(config.sensitivity_para(index,1))*ratio;
                case 5
                    G(config.sensitivity_para(index,1)) = G(config.sensitivity_para(index,1))*ratio;
            end
        end
        tFun = TheoryData(kz, kr, G, d, vhc, w, time, config);
        if config.mode == 'p'
            result(index, :) = (tFun-OriginalFun)./log(ratio);
        else
            result(index, :) = log(tFun./OriginalFun)./log(ratio);
        end
        Sens_Data(:,index+1) = result(index, :)';
    end
    if withInput >= 2
        OutputPath = varargin{2};
    else
        OutputPath = uigetdir('.','Pick a folder to storage the results');
        if isequal(OutputPath,0)
            disp('User pressed cancel')
            fclose all;
            rmdir(temp_path, 's');
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
    mkdir(OutputFolder);
    copyfile(config_file, fullfile(OutputFolder,'configuration.txt'));
    % save the sensitivity picture to the output file
    numPic = length(config.num_line_pic);
    linePic = zeros(numPic,2);
    indexLine = 0;
    for index = 1:1:numPic
        linePic(index, 1) = indexLine + 1;
        indexLine = indexLine + config.num_line_pic(index);
        linePic(index, 2) = indexLine;
    end
    for index = 1:1:numPic
        fig = figure('Position', fPosition);
        hold on
        for indexLine = linePic(index, 1):linePic(index, 2)
            plot(time_ns, result(indexLine, :), 'Color', Colors(indexLine,:),'LineWidth', 2);
        end
        legend(para_name(linePic(index, 1):linePic(index, 2)));
        hold off
        saveas(fig,fullfile(OutputFolder,['Pic_', num2str(index),'.png']),'png');
    end
    Sens_file = fopen(fullfile(OutputFolder, 'Sensitivity_data.txt'),'w+');
    fprintf(Sens_file, '%s\r\n', para_name_s);
    fprintf(Sens_file, form_s, Sens_Data');
    fclose(Sens_file);
end

%% ParameterScan calculation model
if config.ParameterScan == 1
    NumParame = size(config.parameterScan_para, 1);
    Numtime = 200;
    time = linspace(config.tau(1), config.tau(2), Numtime);
    time_ns = time*1e9;
    NumXValue = 50;
    Ratio = 1.01;
    para_name = getLabel_s(config.LayerName, config.parameterScan_para);
    para_name_2 = getLabel_s(config.LayerName, config.parameterScan_para(:,5:6));
    if withInput >= 2
        OutputPath = varargin{2};
    else
        OutputPath = uigetdir('.','Pick a folder to storage the results');
        if isequal(OutputPath,0)
            disp('User pressed cancel')
            fclose all;
            rmdir(temp_path, 's');
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
    mkdir(OutputFolder);
    copyfile(config_file, fullfile(OutputFolder,'configuration.txt'));
    for index_para = 1:1:NumParame
        XValue = linspace(config.parameterScan_para(index_para,3), config.parameterScan_para(index_para,4), NumXValue);
        result = zeros(NumXValue, Numtime);
        para_name_s = 'Time[ns]';
        form_s = '%f\t';
        for indexXValue = 1:1:NumXValue
            para_name_s = [para_name_s,', ',num2str(XValue(indexXValue))];
            form_s = [form_s,'%f\t'];
            kz = config.kz;
            kr = config.kr;
            G = config.G;
            d = config.d;
            vhc = config.vhc;
            w = config.w;
            if config.parameterScan_para(index_para,1) == 0
                switch config.parameterScan_para(index_para,2)
                    case 1
                        w(1) = XValue(indexXValue);
                    case 2
                        w(2) = XValue(indexXValue);
                end
                config_new = config; 
            elseif config.parameterScan_para(index_para,1) == -1
                config_new = config; 
                config_new.f_mod = XValue(indexXValue);
            else
                switch config.parameterScan_para(index_para,2)
                    case 1
                        kz(config.parameterScan_para(index_para,1)) = XValue(indexXValue);
                    case 2
                        kr(config.parameterScan_para(index_para,1)) = XValue(indexXValue);
                    case 3
                        vhc(config.parameterScan_para(index_para,1)) = XValue(indexXValue);
                    case 4
                        d(config.parameterScan_para(index_para,1)) = XValue(indexXValue);
                    case 5
                        G(config.parameterScan_para(index_para,1)) = XValue(indexXValue);
                end
                config_new = config; 
            end
            OriginalFun = TheoryData(kz, kr, G, d, vhc, w, time, config_new);
            if config.parameterScan_para(index_para,5) == 0
                switch config.parameterScan_para(index_para,6)
                    case 1
                        w(1) = w(1)*Ratio;
                    case 2
                        w(2) = w(2)*Ratio;
                end
            elseif config.parameterScan_para(index_para,5) == -1
                config_new.f_mod = config_new.f_mod*Ratio;
            else
                switch config.parameterScan_para(index_para,6)
                    case 1
                        kz(config.parameterScan_para(index_para,5)) = kz(config.parameterScan_para(index_para,5))*Ratio;
                    case 2
                        kr(config.parameterScan_para(index_para,5)) = kr(config.parameterScan_para(index_para,5))*Ratio;
                    case 3
                        vhc(config.parameterScan_para(index_para,5)) = vhc(config.parameterScan_para(index_para,5))*Ratio;
                    case 4
                        d(config.parameterScan_para(index_para,5)) = d(config.parameterScan_para(index_para,5))*Ratio;
                    case 5
                        G(config.parameterScan_para(index_para,5)) = G(config.parameterScan_para(index_para,5))*Ratio;
                end
            end
            tFun = TheoryData(kz, kr, G, d, vhc, w, time, config_new);
            if config.mode == 'p'
                result(indexXValue, :) = (tFun-OriginalFun)./log(Ratio);
            else
                result(indexXValue, :) = log(tFun./OriginalFun)./log(Ratio);
            end
        end
        form_s = [form_s,'\n'];
        [XValue_M, Time_M] = meshgrid(time_ns,XValue);
        fig = figure('Position', fPosition);
        hold on
        contourf(XValue_M,Time_M,result,12,'LineStyle','none');
        colorbar
        xlabel('Time [ns]');
        ylabel(para_name{index_para});
        title(['The sensitivity of ', para_name_2{index_para},' by scanning ', para_name{index_para}]);
        hold off
        saveas(fig,fullfile(OutputFolder,[para_name_2{index_para},'__',para_name{index_para},'.png']),'png');
        Data_file = fopen(fullfile(OutputFolder, [para_name_2{index_para},'__', para_name{index_para},'.txt']),'w+');
        fprintf(Data_file, '%s\r\n', para_name_s);
        data_all = [time_ns' result'];
        fprintf(Data_file, form_s, data_all');
        fclose(Data_file);
    end
end

%% Electron-Phonon Coupling factor
if config.ElectronPhononMode == 1
    length_para = size(config.fit_para,1);
    [para_name, format_f] = getLabel_EP(config.fit_para);
    % all data files in the selected folder will be processed
    if config.folder_mode == 1
        % the folder or file path of the data
        if withInput >= 2
            SourcePath = varargin{2};
        else
            SourcePath = uigetdir('.','Pick a folder where there are data files');
            if isequal(SourcePath,0)
                disp('User pressed cancel')
                fclose all;
                rmdir(temp_path, 's');
                return
            end
        end
        if withInput >= 3
            OutputPath = varargin{3};
        else
            OutputPath = uigetdir('.','Pick a folder to storage the results');
            if isequal(OutputPath,0)
                disp('User pressed cancel')
                fclose all;
                rmdir(temp_path, 's');
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
            fclose all;
            rmdir(temp_path, 's');
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
            Result = TDTR_EP_DataFitting(raw_data, config);
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
            plot(Result.dealed_data.tau,Result.dealed_data.fun,'ko','MarkerSize',8);
            title(['Fitting: ', num2str(index)]);
            hold on
            plot(Result.theory_data.tau,Result.theory_data.fun,'r-','LineWidth',4);
            if config.mode ~= 'p'
                ylim([0,max(Result.dealed_data.fun)+0.5]) 
            end
            legend('Data','Best fit')
            xlabel('time [ps]')
            switch config.mode
                case 'x'
                    ylabel('In-phase signal')
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
                fclose all;
                rmdir(temp_path, 's');
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
                fclose all;
                rmdir(temp_path, 's');
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
        Result = TDTR_EP_DataFitting(raw_data, config);
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
        plot(Result.dealed_data.tau,Result.dealed_data.fun,'ko','MarkerSize',8);
        title(['Fitting: ', num2str(index)]);
        hold on
        plot(Result.theory_data.tau,Result.theory_data.fun,'r-','LineWidth',4);
        if config.mode ~= 'p'
            ylim([0,max(Result.dealed_data.fun)+0.5]) 
        end
        legend('Data','Best fit')
        xlabel('Delay time [ps]')
        switch config.mode
            case 'x'
                ylabel('In-phase signal')
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

%% Mapping
if config.mapping_mode == 1
    length_para = size(config.fit_para,1);
    [para_name, format_f] = getLabel(config.LayerName, config.fit_para);
    % all data files in the selected folder will be processed

    % the folder or file path of the data
    if withInput >= 2
        SourcePath = varargin{2};
    else
        SourcePath = uigetdir('.','Pick a folder where there are data files');
        if isequal(SourcePath,0)
            disp('User pressed cancel')
            fclose all;
            rmdir(temp_path, 's');
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
            fclose all;
            rmdir(temp_path, 's');
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
    filelist = dir([SourcePath filesep() '*.txt']);
    length_filelist = length(filelist);
    if length_filelist < 1
        disp('There is not data file in selected folder.');
        fclose all;
        rmdir(temp_path, 's');
        return
    end
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
    % To do some modification of the origina data
    if config.modification_mode == 1
        [TDTRFileName, TDTRFileDir] = uigetfile('.txt','Pick a TDTR data file, which will be used to modify the Offset data');
        TDTRFilePath = fullfile(TDTRFileDir, TDTRFileName);
        if isequal(TDTRFilePath,0)
            disp('User pressed cancel')
            fclose all;
            rmdir(temp_path, 's');
            return
        end
        raw_data = load(TDTRFilePath);
        [Delta_tau, Delta_phase] = dealTDTRRawData(raw_data, config.ZeroPointMode);
        config.modification_value = [Delta_tau, Delta_phase];
    end
    % To deal the index of the results
    [X_max, Y_max, Index_list] = dealFilelist(filelist);
    Results = zeros(length_filelist,length_para);
    Results_matrix = zeros(X_max, Y_max, length_para);
    for index = 1:1:length_filelist
        disp(['File ', num2str(index), ' : ', filelist(index).name]);
        fprintf(result_file,'%5s %3s %3s %s\r\n','File', num2str(index), ':', filelist(index).name);
        % read raw data
        raw_data = load(fullfile(SourcePath, filelist(index).name));
        % do the fitting and get result
        Result = TDTRDataFitting_Mapping(raw_data, config);
        % save the raw data file to output folder
        raw_data_file = fullfile(raw_data_folder,filelist(index).name);
        copyfile(fullfile(SourcePath, filelist(index).name), raw_data_file);
        % save the dealed data file to output folder
        dealed_data_file = fopen(fullfile(dealed_data_folder, filelist(index).name),'w+');
        fprintf(dealed_data_file, '%f\t%f\r\n', [Result.dealed_data.tau(:)'; Result.dealed_data.fun(:)']);
        fclose(dealed_data_file);
        if index == 1
            Data_length = size(Result.dealed_data.tau(:),1);
            Data_Matrix = zeros(X_max, Y_max, Data_length);
            Time_list = Result.dealed_data.tau(:)';
        end
        Data_Matrix(Index_list(index,1),Index_list(index,2),:) = Result.dealed_data.fun(:);
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
        close(fig);
        % diaplay the result of fitting
        disp(para_name);
        disp(Result.fittingValue);
        disp(['The standard deviation is ', num2str(Result.StdDev)]);
        % save the fitting result to result_log file
        fprintf(result_file,'%s\r\n',para_name);
        fprintf(result_file,format_f,Result.fittingValue');
        fprintf(result_file,'%s\r\n',['The standard deviation is ', num2str(Result.StdDev)]);
        Results(index,:) = Result.fittingValue;
        Results_matrix(Index_list(index,1),Index_list(index,2),:) = Result.fittingValue;
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
    para_name_list = getLabel_s(config.LayerName, config.fit_para);
    % make the image folder
    mapping_folder = fullfile(OutputFolder,'mapping');
    mkdir(mapping_folder);
    mapping_parameter_folder = fullfile(mapping_folder,'parameters');
    mkdir(mapping_parameter_folder);
    for index = 1:1:length_para
        para_data = Results_matrix(:,:,index);
        length_x = size(para_data,2);
        length_y = size(para_data,1);
        para_filepath = fullfile(mapping_parameter_folder,[para_name_list{index} '.txt']);
        save(para_filepath,'para_data','-ascii');
        fig = figure('Position', fPosition);
        imagesc([0,length_x]*config.interval(1), [0,length_y]*config.interval(2), para_data);
        axis equal;
        colorbar
        title(para_name_list{index});
        saveas(fig,fullfile(mapping_parameter_folder,[para_name_list{index},'.png']),'png');
    end
    mapping_data_folder = fullfile(mapping_folder, 'data');
    mkdir(mapping_data_folder);
    for index = 1:1:Data_length
        data = Data_Matrix(:,:,index);
        length_x = size(data,2);
        length_y = size(data,1);
        switch config.mode
            case 'r'
                data_type = 'Ratio';
            case 'p'
                data_type = 'Phase';
            case 'a'
                data_type = 'Amplitude';
        end
        data_filepath = fullfile(mapping_data_folder,['time-' num2str(Time_list(index)) 'ns(' data_type ').txt']);
        save(data_filepath,'data','-ascii');
        fig = figure('Position', fPosition);
        imagesc([0,length_x]*config.interval(1), [0,length_y]*config.interval(2), data);
        axis equal;
        colorbar
        title(['time:' num2str(Time_list(index)) 'ns(' data_type ')']);
        saveas(fig,fullfile(mapping_data_folder,['time-' num2str(Time_list(index)) 'ns(' data_type ').png']),'png');
        close(fig);
    end
end

%% close file and delete the temp folder
fclose all;
rmdir(temp_path, 's');