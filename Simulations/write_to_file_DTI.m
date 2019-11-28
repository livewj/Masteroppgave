function write_to_file_DTI(foldername, filename, FA_app, MD_app, RD_app, AD_app)

fileID = fopen([foldername filename],'w');

%fprintf(fileID, '%2s \n', 'FA_app');
for i=1:length(FA_app)
    fprintf(fileID, '%18.5f', FA_app(i,:));
    fprintf(fileID, '\n');
end

%fprintf(fileID, '%2s \n', 'MD_app');
for i=1:length(MD_app)
    fprintf(fileID, '%18.5f', MD_app(i,:));
    fprintf(fileID, '\n');
end

%fprintf(fileID, '%2s \n', 'RD_app');
for i=1:length(RD_app)
    fprintf(fileID, '%18.5f', RD_app(i,:));
    fprintf(fileID, '\n');
end

%fprintf(fileID, '%2s \n', 'AD_app');
for i=1:length(AD_app)
    fprintf(fileID, '%18.5f', AD_app(i,:));
    fprintf(fileID, '\n');
end

fclose(fileID);


end