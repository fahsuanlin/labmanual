close all; clear all;

root_path='/space/maki8/1/fhlin/seeg_language';
subject={
's041';
's046';
's050';
's052';
's054';
};

fstem={
'seeg_wb_mne_091019_l_mne';
'seeg_wb_mne_091019_n_mne';
};

target_subject='fsaverage';
pdir=pwd;

for subj_idx=1:length(subject)
	cd(sprintf('%s/%s/analysis',root_path,subject{subj_idx}));


        %morphing STC....
	for f_idx=1:length(fstem)
		for hemi_idx=1:2
			switch hemi_idx
			case 1
                                hemi='lh';
                        case 2
                                hemi='rh';
                        end;
                        fn_in=sprintf('%s-%s.stc',fstem{f_idx},hemi);
                        fn_out=sprintf('%s_2_%s_%s',subject{subj_idx},target_subject,fstem{f_idx});
                        cmd=sprintf('!mne_make_movie --subject %s --stcin %s --morph %s --stc %s --%s --smooth 5', subject{subj_idx}, fn_in, target_subject, fn_out, hemi);
                        eval(cmd);
                end;
        end;
        %fn=sprintf('%s_%s-vol.stc',output_stem,erf_all(trig_idx).trig_str);
        %fprintf('\tsaving [%s]...\n',fn);
        %inverse_write_stc(X_dspm,[0:size(X_dspm,1)-1],min(erf_all(trig_idx).timeVec),t0,fn);

end;
cd(pdir);
fprintf('Done!\n');

