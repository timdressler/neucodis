% pp_subj_sanity_new.m
%
% Perform subject sanity check.
% Using preprocessed data.
% Preprocessing done with pp_erp_pre_proc.m & pp_coherence_pre_proc.m.
% Check if subjects had enough valid trials for coherence & ERP analysis.
% Stores subjects with enough valid trials in separate folder.
% Subjects were removed from analysis if they had less than 50 trials in
% *either* preprocessing script. Therefore, the subjects analyzed are the
% same for both ERP and coherence analysis.
%
% Tim Dressler, 18.11.2024