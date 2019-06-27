%> @file  SCPLearnFromMatFiles.m
%> @brief Function that runs SCPLearn on nifti time-series data 
%======================================================================
%> @brief It takes as input subject level mat files and runs SCPLearn
%> It returns 
%> (1) all ROI timeseries in a mat file
%> (2) results of SCPLearn in a mat file
%> For details see the following paper:
%> http://www.sciencedirect.com/science/article/pii/S1053811914008003
%>
%> @param matlist Text file containing list of mat files
%> @param K Number of SCPs
%> @param pruning SCPs with inner-product overlap > pruning are discarded
%> @param lambda Sparsity level
%> @param outprefix prefix for all output files (will be overwritten if exists)
%> @param levels {0,1}, if 1 will run heirarchical SCP decomposition
%> @param verbose Integer value if 1 verbose messages will be output
%>
%> @b Author: 
%> Harini Eavani
%>
%> @b Link: 
%> https://www.cbica.upenn.edu/sbia/software/
%> 
%> @b Contact: 
%> sbia-software@uphs.upenn.edu
%======================================================================
function [] = SCPLearnFromMatFiles(matlist,K,pruning,lambda,outprefix,levels,verbose,test)



matfile=[outprefix,'_ts.mat'];
scpmatfile=[outprefix,'_SCPs.mat'];
csvfile=[outprefix,'_SCP_Coeffs.csv'];
basisfile=[outprefix,'_SCP_basis.csv'];

%> open datafile
fid = fopen(matlist,'r');
if(fid<0)
    %> if cant open return with error msg
    fprintf('Count not open %s\n',matlist)
    return
end

%> read first line
tline = fgetl(fid);
count=1;
sample_weights=[];
ids={};
ts={};
while ischar(tline) %> while read line is a char
    [filename,fields]  = strtok(tline,',');
    [fields,~]  = strtok(fields,',');
    
    %> load nifti file
    if(verbose)
        fprintf('Loading %s\n',filename)
    end
    try
        data_var = load(filename);
    catch err
        fprintf('Could not load %s, are inputs nifti?\n',filename)
        tline = fgetl(fid);
        count = count + 1;
        continue
    end
    names = fieldnames(data_var);
    ts{count} = data_var.(names{1});

	% normalize
	ts{count} = bsxfun(@minus,ts{count},mean(ts{count},2));
    ts{count} = bsxfun(@times,ts{count},1./std(ts{count},[],2));
    
    sample_weights=[sample_weights;str2double(fields)];
    ids{count}=filename;
    tline = fgetl(fid);
    count = count + 1;
    
end
fclose(fid); %> close txt file

%> save time-series in mat file
if(~any(isnan(sample_weights)))
    sample_weights=reshape(sample_weights,length(sample_weights),1);
    save(matfile,'ts','sample_weights')
else
    save(matfile,'ts')
end

%> Run SCPLearn
vars = whos('-file',matfile);
if(ismember('sample_weights', {vars.name}))
    fprintf('Found sample weights as input\n')
    return
else
    if (test == '0')
        SCPLearn(matfile, K, lambda, scpmatfile,verbose,pruning,levels)
    else
        SCPTest(matfile,test,scpmatfile,verbose)
    end
end

load(scpmatfile)
%> write csv file
fp=fopen(csvfile,'w');
header_string=['ID'];
for kk=1:size(C,1)
    header_string=[header_string,',SCP_',num2str(kk)];
end
fprintf(fp,'%s\n',header_string);

for ii=1:length(ids)
    temp = strsplit(ids{ii},'/');
    string=[temp{end}];
    for kk=1:size(C,1)
        string=[string,',',num2str(C(kk,ii))];
    end
    fprintf(fp,'%s\n',string);
end
fclose(fp);    

%> write basis file
if (exist('B','var'))
    fp=fopen(basisfile,'w');
    header_string=['ROI_name'];
    for kk=1:size(B,2)
        header_string=[header_string,',SCP_',num2str(kk)];
    end
    fprintf(fp,'%s\n',header_string);

    for ii=1:size(B,1)
        string=[num2str(ii)];
        for kk=1:size(B,2)
            string=[string,',',num2str(B(ii,kk))];
        end
        fprintf(fp,'%s\n',string);
    end
    fclose(fp); 
end    
end
