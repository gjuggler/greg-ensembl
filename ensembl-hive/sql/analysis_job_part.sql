------------------------------------------------------------------------------------
--
-- Table structure for table 'analysis_job'
--
-- overview:
--   The analysis_job is the heart of this sytem.  It is the kiosk or blackboard
--   where workers find things to do and then post work for other works to do.
--   The job_claim is a UUID set with an UPDATE LIMIT by worker as they fight
--   over the work.  These jobs are created prior to work being done, are claimed
--   by workers, are updated as the work is done, with a final update on completion.
--
-- semantics:
--   analysis_job_id         - autoincrement id
--   prev_analysis_job_id    - previous analysis_job which created this one (and passed input_id)
--   analysis_id             - the analysis_id needed to accomplish this job.
--   input_id                - input data passed into Analysis:RunnableDB to control the work
--   job_claim               - UUID set by workers as the fight over jobs
--   hive_id                 - link to hive table to define which worker claimed this job
--   status                  - state the job is in
--   retry_count             - number times job had to be reset when worker failed to run it
--   completed               - timestamp when job was completed
--   branch_code             - switch-like branching control, default=1 (ie true)

DROP TABLE IF EXISTS analysis_job_partit;
-- CREATE TABLE analysis_job_partit (
--   analysis_job_id           int(10) NOT NULL auto_increment,
--   prev_analysis_job_id      int(10) NOT NULL,  #analysis_job which created this from rules
--   analysis_id               int(10) NOT NULL,
--   input_id                  char(255) not null,
--   job_claim                 char(40) NOT NULL default '', #UUID
--   hive_id                   int(10) NOT NULL,
--   status                    enum('READY','BLOCKED','CLAIMED','GET_INPUT','RUN','WRITE_OUTPUT','DONE','FAILED') DEFAULT 'READY' NOT NULL,
--   retry_count               int(10) default 0 not NULL,
--   completed                 datetime NOT NULL,
--   branch_code               int(10) default 1 NOT NULL,
--   runtime_msec              int(10) default 0 NOT NULL, 
--   query_count               int(10) default 0 NOT NULL, 

--   PRIMARY KEY                 (analysis_id, analysis_job_id),
--   UNIQUE KEY                  (analysis_job_id),
-- #  UNIQUE KEY input_id_analysis (input_id, analysis_id),
-- #  UNIQUE KEY input_id_analysis (analysis_id, input_id),
--   INDEX claim_analysis_status  (job_claim, analysis_id, status),
--   INDEX analysis_status        (analysis_id, status),
--   INDEX hive_id                (hive_id)
-- ) ENGINE=InnoDB MAX_ROWS = 100000000 AVG_ROW_LENGTH = 390 PARTITION BY KEY(analysis_id) PARTITIONS 100;


CREATE TABLE analysis_job_partit (
  analysis_job_id           int(10) NOT NULL auto_increment,
  prev_analysis_job_id      int(10) NOT NULL,  #analysis_job which created this from rules
  analysis_id               int(10) NOT NULL,
  input_id                  char(255) not null,
  job_claim                 char(40) NOT NULL default '', #UUID
  hive_id                   int(10) NOT NULL,
  status                    enum('READY','BLOCKED','CLAIMED','GET_INPUT','RUN','WRITE_OUTPUT','DONE','FAILED') DEFAULT 'READY' NOT NULL,
  retry_count               int(10) default 0 not NULL,
  completed                 datetime NOT NULL,
  branch_code               int(10) default 1 NOT NULL,
  runtime_msec              int(10) default 0 NOT NULL, 
  query_count               int(10) default 0 NOT NULL, 

  PRIMARY KEY                  (analysis_job_id),
  # UNIQUE KEY input_id_analysis (input_id, analysis_id),
  INDEX claim_analysis_status  (job_claim, analysis_id, status),
  INDEX analysis_status        (analysis_id, status),
  INDEX hive_id                (hive_id)
) ENGINE=InnoDB MAX_ROWS = 100000000 AVG_ROW_LENGTH = 390 PARTITION BY KEY(analysis_id) PARTITIONS 100;

-- CREATE TABLE analysis_job (
--   analysis_job_id           int(10) NOT NULL auto_increment,
--   prev_analysis_job_id      int(10) NOT NULL,  #analysis_job which created this from rules
--   analysis_id               int(10) NOT NULL,
--   input_id                  char(255) not null,
--   job_claim                 char(40) NOT NULL default '', #UUID
--   hive_id                   int(10) NOT NULL,
--   status                    enum('READY','BLOCKED','CLAIMED','GET_INPUT','RUN','WRITE_OUTPUT','DONE','FAILED') DEFAULT 'READY' NOT NULL,
--   retry_count               int(10) default 0 not NULL,
--   completed                 datetime NOT NULL,
--   branch_code               int(10) default 1 NOT NULL,
--   runtime_msec              int(10) default 0 NOT NULL, 
--   query_count               int(10) default 0 NOT NULL, 

--   PRIMARY KEY                  (analysis_job_id),
--   UNIQUE KEY input_id_analysis (input_id, analysis_id),
--   INDEX claim_analysis_status  (job_claim, analysis_id, status),
--   INDEX analysis_status        (analysis_id, status),
--   INDEX hive_id                (hive_id)
-- ) ENGINE=InnoDB MAX_ROWS = 100000000 AVG_ROW_LENGTH = 390 ; 

