DELETE FROM analysis_job WHERE analysis_id=4;
UPDATE analysis_stats SET hive_capacity=30 WHERE analysis_id=4;

set @primates="human,chimpanzee,gorilla,orangutan,macaque,marmoset,mouse";

insert into analysis_job (analysis_id,status,input_id) VALUES

  # Remember: by default we have model=2 (foreground-model) and foreground_species => '9593,9606'
#  (4,"READY",CONCAT("{data_id=>'999',job_role=>'fan_jobs',keep_species=>'",@primates,"'}"))

  (4,"READY",CONCAT("{data_id=>'101',target_alignment_length=>'100000',branch_model=>1,keep_species=>'",@primates,"'}")),
  (4,"READY",CONCAT("{data_id=>'102',target_alignment_length=>'500000',branch_model=>1,keep_species=>'",@primates,"'}")),
  (4,"READY",CONCAT("{data_id=>'103',target_alignment_length=>'1000000',branch_model=>1,keep_species=>'",@primates,"'}")),
  (4,"READY",CONCAT("{data_id=>'104',target_alignment_length=>'5000000',branch_model=>1,keep_species=>'",@primates,"'}")),
  (4,"READY",CONCAT("{data_id=>'105',target_alignment_length=>'10000000',branch_model=>1,keep_species=>'",@primates,"'}"))

;