dataset_path = '\\zserver.curie.fr\umr3664\equipe_dostatni\g_fernandes\RAW\190509_GF1.2\table_summary';
dataset_file = 'th35_21'; % No extension here

change_fix = 1;
backup = 1;

offset=120;
try
    datamat=dlmread(fullfile(dataset_path,[dataset_file '.txt']),' ');
catch
    datamat=dlmread(fullfile(dataset_path,[dataset_file '.txt']),',');
end
dlmwrite(fullfile(dataset_path,[dataset_file '_backup.txt']),datamat);
datamat(:,16)=datamat(:,16)-120;
datamat(:,17)=datamat(:,17)+120;
dlmwrite(fullfile(dataset_path,[dataset_file '.txt']),datamat,' ');
display('Original file corrected with new shiftL and shiftR');
if change_fix
    if exist(fullfile(dataset_path,[dataset_file '_fixed.txt']),'file')
        % Read
        try
            datamat=dlmread(fullfile(dataset_path,[dataset_file '_fixed.txt']),' ');
        catch
            datamat=dlmread(fullfile(dataset_path,[dataset_file '_fixed.txt']),',');
        end
        % Backup
        dlmwrite(fullfile(dataset_path,[dataset_file '_fixed_backup.txt']),datamat);
        % Fix
        datamat(:,16)=datamat(:,16)-120;
        datamat(:,17)=datamat(:,17)+120;
        % Overwrite        
        dlmwrite(fullfile(dataset_path,[dataset_file '_fixed.txt']),datamat);
        display('Fixed file corrected with new shiftL and shiftR');
    else
        display('No fixed file found');
    end
end
