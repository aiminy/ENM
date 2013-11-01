function batch_calculate_random_rmsip(n,pdb_file)

 fid=fopen('anm_random_pca_rmsip_5853.txt','a');
 
 for i=1:n
  [rmsip_pca_anm,rmsip_pca_random,anm_random]=calculate_random_anm_pc_rmsip(pdb_file);
  
  pc_random(i)=rmsip_pca_random;
  anm_ran(i)=anm_random;
  
  
 end

fprintf(fid,'%8s %6.4f %6.4f %6.4f \n',char(pdb_file),rmsip_pca_anm,mean(pc_random),mean(anm_ran));
fclose(fid);
