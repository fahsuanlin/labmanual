close all; clear all;

subject='s015';

eval(sprintf('!om_assemble -HeadMat %s_d10_bem.geom %s_d10_bem.cond %s_d10.hm',subject,subject,subject));
eval(sprintf('!om_assemble -DSM %s_d10_bem.geom %s_d10_bem.cond %s_d10-lh.dip %s_d10-lh.dsm brain',subject,subject,subject,subject));
eval(sprintf('!om_assemble -DSM %s_d10_bem.geom %s_d10_bem.cond %s_d10-rh.dip %s_d10-rh.dsm brain',subject,subject,subject,subject));
eval(sprintf('!om_assemble -h2em %s_d10_bem.geom %s_d10_bem.cond %s_d10.eegsensors %s_d10.h2em',subject,subject,subject,subject));
eval(sprintf('!om_minverser %s_d10.hm %s_d10.hm_inv',subject,subject));
eval(sprintf('!om_gain -EEG %s_d10.hm_inv %s_d10-lh.dsm %s_d10.h2em eeg_lh_gain.mat',subject,subject,subject));
eval(sprintf('!om_gain -EEG %s_d10.hm_inv %s_d10-rh.dsm %s_d10.h2em eeg_rh_gain.mat',subject,subject,subject));

