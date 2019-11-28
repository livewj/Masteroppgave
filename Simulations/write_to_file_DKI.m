function write_to_file_DKI(foldername, filename, fa_app, md_app, rd_app, ad_app, mk_app, rk_app, ak_app)

fileID = fopen([foldername filename],'w');

%fprintf(fileID, '%2s \n', 'fa_app');
for i=1:length(fa_app)
    fprintf(fileID, '%18.5f', fa_app(i,:));
    fprintf(fileID, '\n');
end

%fprintf(fileID, '%2s \n', 'md_app');
for i=1:length(md_app)
    fprintf(fileID, '%18.5f', md_app(i,:));
    fprintf(fileID, '\n');
end

%fprintf(fileID, '%2s \n', 'rd_app');
for i=1:length(rd_app)
    fprintf(fileID, '%18.5f', rd_app(i,:));
    fprintf(fileID, '\n');
end

%fprintf(fileID, '%2s \n', 'ad_app');
for i=1:length(ad_app)
    fprintf(fileID, '%18.5f', ad_app(i,:));
    fprintf(fileID, '\n');
end

%fprintf(fileID, '%2s \n', 'mk_app');
for i=1:length(mk_app)
    fprintf(fileID, '%18.5f', mk_app(i,:));
    fprintf(fileID, '\n');
end

%fprintf(fileID, '%2s \n', 'rk_app');
for i=1:length(rk_app)
    fprintf(fileID, '%18.5f', rk_app(i,:));
    fprintf(fileID, '\n');
end

%fprintf(fileID, '%2s \n', 'ak_app');
for i=1:length(ak_app)
    fprintf(fileID, '%18.5f', ak_app(i,:));
    fprintf(fileID, '\n');
end

fclose(fileID);


end