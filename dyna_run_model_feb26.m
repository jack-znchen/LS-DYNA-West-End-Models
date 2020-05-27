close all
clear

%read building summary spreadsheet
bldg= readmatrix('bldg_table_removed small area.xlsx');

%matrix to store 
summary= zeros(size(bldg,1),2);

for bldg_num= 1:size(bldg,1)
    %extract information needed to create .key file
    bldg_id= bldg(bldg_num,21);
    summary(bldg_num,2)= bldg(bldg_num,4);
    
    %run lsdyna using command prompt
    directory1= 'cd\Users\billc\Downloads\Matlab\Models_Feb26\bldg';
    directory2= '_file & C:\LSDYNA\ls-dyna_smp_d_R920_winx64_ifort131.exe i=C:\Users\billc\Downloads\Matlab\Models_Feb26\bldg';
    directory3= '_file\Bldg';
    
    command= sprintf('%s%d%s%d%s%d_Keyword_File.key',directory1,bldg_id,directory2,bldg_id,directory3,bldg_id);
    status = system(command);
    
    %turn eigout file into a text file
    directory4= 'rename  C:\Users\billc\Downloads\Matlab\Models_Feb26\bldg';
    directory5= '_file\eigout   eigout.txt';
    
    system(sprintf('%s%d%s',directory4,bldg_id,directory5));
    
    %read eigout file
    directory6= 'Models_Feb26\bldg';
    directory7= '_file\eigout.txt';
    dyna_file =fopen(sprintf('%s%d%s',directory6,bldg_id,directory7),'r');
    
    i = 1;
    line = fgetl(dyna_file);
    read_file{i,1} = line;

    while ischar(line)
        i = i+1;
        line = fgetl(dyna_file);
        read_file{i,1} = line;
    end

    fclose(dyna_file);
    
    %read and store fundamental period
    A = read_file{12};
    B= str2num(A);
    summary(bldg_num,1)= B(end);
end

%plot results
scatter(summary(:,1),summary(:,2))

xlim([0 6])
ylim([0 35])

xlabel('Period (s)');
ylabel('Storys Above Grade');