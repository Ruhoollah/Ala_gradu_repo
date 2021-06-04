
%%  MAIN parameters
% set the input directory where your data is stored
data_dir = 'C:\Users\Ala\Rdata';
% result_dir ='C:\Users\Ala\Alaslab';
% specify the file type of your data
data_type = '*.set';
% use sbj_filt to select all (or a subset) of available recordings
sbj_filt = 1 ; %setdiff(can1:12, [3 7]);
% use ctapID to uniquely name the base folder of the output directory tree
ctapID = 'Ala_blink_result';
% use keyword 'all' to select all stepSets, or use some index
set_select = 'All';
% set the electrode for which to calculate and plot ERPs after preprocessing
erploc = 'A1';

% Runtime options for CTAP:
STOP_ON_ERROR = true;
OVERWRITE_OLD_RESULTS = true;

%% Create the CONFIGURATION struct
% First, define step sets & their parameters: sbf_cfg() is written by the USER
[Cfg, ctap_args] = sbf_cfg('C:\Users\Ala', ctapID);

% Select step sets to process
Cfg.pipe.runSets = set_select;

% Next, create measurement config (MC) based on folder, & select subject subset
Cfg = get_meas_cfg_MC(Cfg, data_dir, 'eeg_ext', data_type, 'sbj_filt', sbj_filt);

% Assign arguments to the selected functions, perform various checks
Cfg = ctap_auto_config(Cfg, ctap_args);

%% Run the pipe
tic
CTAP_pipeline_looper(Cfg,...
                    'debug', STOP_ON_ERROR,...
                    'overwrite', OVERWRITE_OLD_RESULTS);
toc
%cleanup the global workspace
clear STOP_ON_ERROR OVERWRITE_OLD_RESULTS ctap_args

%% Subfunctions
% Pipe definition
function [Cfg, out] = sbf_cfg(project_root_folder, ID)


%% Define important directories and files
% Analysis ID
Cfg.id = ID;
% Directory where to locate project - in this case, just the same as input dir
Cfg.env.paths.projectRoot = project_root_folder;
% CTAP root dir named for the ID
Cfg.env.paths.ctapRoot = fullfile(Cfg.env.paths.projectRoot, Cfg.id);
% CTAP output goes into analysisRoot dir, here can be same as CTAP root
Cfg.env.paths.analysisRoot = Cfg.env.paths.ctapRoot;
%Cfg.eeg.chanlocs = fullfile('C:\','Users','Ala','EEGdata','Sub_2','chanlocs_sub3','sub-003_FaceRecog_electrodes.tsv');

%% Define other important stuff
Cfg.eeg.reference = {'average'};

% NOTE! EOG channel specification for artifact detection purposes. The
% HeadIT dataset did not use a traditional EOG montage with extra channels
% above and below the eye, as there were EEG channels mounted all the
% way down the forehead. Thus, EXG1-4 were mounted below and beside the
% eyes. We use H24 over EXG2 to approximate a bipolar lead for vertical EOG
Cfg.eeg.heogChannelNames = {'EXG3' 'EXG4'};
Cfg.eeg.veogChannelNames = {'EXG5' 'EXG6'};
%% Configure analysis pipe

%% Define pipeline 
clear('stepSet');

i = 1; %stepSet 1
stepSet(i).funH = { @CTAP_load_data,... 
                    @CTAP_blink2event,...
                    @CTAP_reject_data};
                
stepSet(i).id = [num2str(i) '_load'];

out.fir_filter = struct(...
   'locutoff', 1);


i = i+1;  %stepSet 2
stepSet(i).funH = { @CTAP_run_ica,...
                    @CTAP_peek_data };
stepSet(i).id = [num2str(i) '_ICA'];
stepSet(i).srcID = '';

out.run_ica = struct(...
    'method', 'fastica',...
    'overwrite', true);
out.run_ica.channelTypes = [1 33 52:111 114:119 131 132];

%i = i+1;  %stepSet 3
%stepSet(i).funH = { @CTAP_detect_bad_comps,... %ADJUST for horizontal eye moves
 %                   @CTAP_reject_data,...
  %                  };
%stepSet(i).id = [num2str(i) '_artifact_correction'];

%ut.detect_bad_comps = struct(...
%    'method', {'adjust' 'blink_template'},...
%    'adjustarg', {'horiz' ''});


i = i+1;  %stepSet 3

stepSet(i).funH = { @CTAP_detect_bad_comps,... %detect blink related ICs
                    @CTAP_filter_blink_ica,...@CTAP_detect_bad_comps,...
                    @CTAP_peek_data};

stepSet(i).id = [num2str(i) '_IC_CORRECTION'];

out.detect_bad_comps = struct(...
    'method', 'recu_blink_tmpl');
out.detect_bad_comps.test_pc = 30; 


out.peek_data = struct(...
    'secs', [10 30],... %start few seconds after data starts
    'peekStats', false,... %get statistics for each peek!
    'overwrite', true,...
    'plotAllPeeks', false,...
    'savePeekData', true,...
    'savePeekICA', true);


%% Store to Cfgg
Cfg.pipe.stepSets = stepSet; % return all step sets inside Cfg struct


end

