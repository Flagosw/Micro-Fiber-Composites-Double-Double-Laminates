function mfc_write_json(file,value)
%MFC_WRITE_JSON Write a UTF-8 JSON evidence record without workspace effects.
fid=fopen(file,'w','n','UTF-8');assert(fid>=0,'MFC:Export','Cannot open output.');
closer=onCleanup(@()fclose(fid));fprintf(fid,'%s\n',jsonencode(value,PrettyPrint=true));
end
