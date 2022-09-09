classdef FileManager < handle
% first parameter: the path where data is stored
% second parameter: the reconstruction project name, a new name will result
% re-process of all the data, while an old name will check which file has
% been processed    
    
    properties
        input_dir
        all_data_file
        output_dir 
        cache_path
        computer_name
        
        txt_dir
    end
    
    methods
        function self = FileManager(input_dir,output_dir,synchronization_path,mask)
            % first parameter: the path where data is stored
            % second parameter: the reconstruction project name, a new name will result
            % re-process of all the data, while an old name will check which file has
            % been processed    

            if nargin <= 3
                mask = '*';
            else
                mask = ['*' mask '*'];
            end
            
            
            
            if ~exist(fullfile(synchronization_path,output_dir,'FileManagerCache'), 'dir')
               mkdir(fullfile(synchronization_path,output_dir,'FileManagerCache'))
            end           
            
            self.cache_path = fullfile(synchronization_path,output_dir,'FileManagerCache');
% 
%             if isfile(fullfile(self.cache_path,'ProcessedDataList.mat'))
%                 processed_data_file = importdata(fullfile(self.cache_path,'ProcessedDataList.mat'));
%             else
%                 processed_data_file = [];
%             end            
%             
            self.input_dir = input_dir;

            if ~exist(fullfile(synchronization_path,output_dir,'Process_Output'), 'dir')
               mkdir(fullfile(synchronization_path,output_dir,'Process_Output'))
            end         
            self.txt_dir =fullfile(synchronization_path,output_dir);
            self.output_dir = fullfile(synchronization_path,output_dir,'Process_Output');


            filenameN = [];
            if iscell(input_dir)
                for dir_index = 1:length(input_dir)
                    Dataset = dir(fullfile(input_dir{dir_index},'**',[ mask '.*']));

                    for file_index =1:length(Dataset)
                        if Dataset(file_index).bytes>100E+6
                            filenameN = cat(1,Dataset(file_index),filenameN);
                        end
                    end   
                end
            else 
                Dataset = dir(fullfile(input_dir,'**','*.*'));
                for file_index =1:length(Dataset)
                    if Dataset(file_index).bytes>100E+6
                        filenameN = cat(1,Dataset(file_index),filenameN);
                    end
                end                
            end

            self.all_data_file = filenameN;
            
            [~, sys_name] = system('hostname');
            sys_name(isletter(sys_name)==0)=[];
            self.computer_name = sys_name;
            
            
        end

        function mask_file(self,mask)
            filenameN = [];
            if iscell(self.input_dir)
                for dir_index = 1:length(self.input_dir)
                    Dataset = dir(fullfile(self.input_dir{dir_index},'**',['*' mask '.*']));

                    for file_index =1:length(Dataset)
                        if Dataset(file_index).bytes>100E+6
                            filenameN = cat(1,Dataset(file_index),filenameN);
                        end
                    end   
                end
            else 
                Dataset = dir(fullfile(self.input_dir,'**','*.*'));
                for file_index =1:length(Dataset)
                    if Dataset(file_index).bytes>100E+6
                        filenameN = cat(1,Dataset(file_index),filenameN);
                    end
                end                
            end

            self.all_data_file = filenameN; 
            

        end
        
        function select_file(self,RetinaList)
            select_data_file = [];
            for ind = 1:RetinaList.get_list_length()
                fileID = RetinaList.get_id(ind-1);
                for adindex = 1:length(self.all_data_file)
                    if contains(fileID,self.all_data_file(adindex).name)
                        select_data_file = cat(1,self.all_data_file(adindex),select_data_file);
                    end
                end  
            end
            for ind = 1:length(select_data_file)
                select_data_file(ind).name
            end
            self.all_data_file = select_data_file;
        end
        
        function log_txt(self,text)
            logID = fopen(fullfile(self.txt_dir,'processing_log.txt'),'a');
            fprintf(logID,[self.computer_name '   ']);
            
            fprintf(logID,datestr(now));
            fprintf(logID,['     ' text '\n']);
            fclose(logID);
        
        end
        
        
        function file_list = get_unprocessed_data(self)
            
            [~, sys_name] = system('hostname');
            sys_name(isletter(sys_name)==0)=[];
            sys_name = [sys_name '.mat'];
            if isfile(fullfile(self.cache_path,sys_name))
                file_list = importdata(fullfile(self.cache_path,sys_name));
            else
                file_list = [];    
            end
            
            if isfile(fullfile(self.cache_path,'ProcessedDataList.mat'))
                processed_data_file = importdata(fullfile(self.cache_path,'ProcessedDataList.mat'));
            else
                processed_data_file = [];
            end 

            for findex = 1:length(self.all_data_file)
                find_it = 0;
                for pfindex = 1:length(processed_data_file)
                    if self.all_data_file(findex).datenum == processed_data_file(pfindex).datenum
                        find_it = 1;
                    end
                end
                if find_it == 0
                    file_list = cat(1,file_list,self.all_data_file(findex));
                end
            end
            %% log
            logID = fopen(fullfile(self.txt_dir,'processing_log.txt'),'a');
            fprintf(logID,[self.computer_name '   ']);
            
            fprintf(logID,datestr(now));
            fprintf(logID,'    %d datasets to process\n',length(file_list));
            fclose(logID);
 
        end
        
        function file_processing(self,file_pointer)
            %% log
            logID = fopen(fullfile(self.txt_dir,'processing_log.txt'),'a');
            fprintf(logID,[self.computer_name '   ']);
            fprintf(logID,datestr(now));
            fprintf(logID,'   begin to process %s_%s \n',file_pointer.folder, file_pointer.name);
            fclose(logID);
            %% add to finish list
            if isfile(fullfile(self.cache_path,'ProcessedDataList.mat'))
                processed_data_file = importdata(fullfile(self.cache_path,'ProcessedDataList.mat'));
            else
                processed_data_file = [];
            end 
            processed_data_file = cat(1,processed_data_file,file_pointer);
            save(fullfile(self.cache_path,'ProcessedDataList.mat'),'processed_data_file');
            
            %% add to ongoing list
            [~, sys_name] = system('hostname');
            sys_name(isletter(sys_name)==0)=[];
            sys_name = [sys_name '.mat'];
            save(fullfile(self.cache_path,sys_name),'file_pointer');
            
            
        end
        
        function set = select_dataset(self,pattern)
            set = [];
            for findex = 1:length(self.all_data_file)
                if contains(self.all_data_file(findex).name,pattern)
                    set = cat(1,set,self.all_data_file(findex));
                end
            end
        end
        
        function file_processed(self,file_pointer)
            %% log

            logID = fopen(fullfile(self.txt_dir,'processing_log.txt'),'a');
            fprintf(logID,[self.computer_name '   ']);
            fprintf(logID,datestr(now));
            fprintf(logID,'   finish %s_%s \n',file_pointer.folder, file_pointer.name);
            fclose(logID);
            
            %% delete ongoing
            [~, sys_name] = system('hostname');
            sys_name(isletter(sys_name)==0)=[];
            sys_name = [sys_name '.mat'];
            delete(fullfile(self.cache_path,sys_name));
            

            
        end            
            
        
        function o = get_output_dir(self)
            o = self.output_dir;
        end
        
        
        function export_file_excel(self,path)
            T = struct2table(self.all_data_file);
            sortedT = sortrows(T,'datenum','ascend'); 
            fl = sortedT.('folder');
            fn = sortedT.('name');

            
            sortedT.IRB = (cellfun(@get_IRB,fl,fn,'UniformOutput',0));
            sortedT.Time = sortedT.date;
            sortedT.Dataset = sortedT.name;
            sortedT.Comment(:) = {' '};%cell(height(sortedT),1);
            sortedT.Quality(:) = {' '};%cell(height(sortedT),1);
            sortedT.Diagnosis(:) = {' '}; %cell(height(sortedT),1);
            sortedT.Protocal(:) = {' '};%cell(height(sortedT),1);
            sortedT = removevars(sortedT,{'name','folder','date','bytes','isdir','datenum'});
            writetable(sortedT,[path '/datalist.xlsx'],'Sheet',1);
            
            %time_sort_index = sort(self.all_data_file)
        end

    end
end


function s = get_IRB(p,n)

    if n(end) == 'x' || n(end-1) == 'x'
        [aa,~,~] = fileparts(p);
        [~,aaa,~] = fileparts(aa);
        s=aaa(1:end-8);
    else
        [~,aaa,~] = fileparts(p);
        s=aaa(1:end-8);
    end

end

