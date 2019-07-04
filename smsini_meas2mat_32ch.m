close all; clear all;

addpath(genpath('/autofs/space/maki6_001/users/eva/toolbox/tool_mb_recon/'));

pdir=pwd;
subfolder_for_meas='32ch';
run_recon_mb_vbvd('/space/maki7/users/fhlin/smsini_array_nccu/063019/meas/32ch_TWIX/',sprintf('%s/%s/',pdir,subfolder_for_meas));

