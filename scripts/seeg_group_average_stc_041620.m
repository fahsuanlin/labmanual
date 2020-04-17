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

output_stem={
'average_041620_seeg_wb_mne_091019_l_z';
'average_041620_seeg_wb_mne_091019_n_z';
};

pdir=pwd;

for f_idx=1:length(fstem)
        for hemi_idx=1:2
                switch hemi_idx
                case 1  
                        hemi='lh';
                case 2
                        hemi='rh';
                end;

		stc=[];

		for subj_idx=1:length(subject)
			cd(sprintf('%s/%s/analysis',root_path,subject{subj_idx}));


                        fn=sprintf('%s_2_%s_%s-%s.stc',subject{subj_idx},target_subject,fstem{f_idx},hemi);
                	[stc(:,:,subj_idx),v,a,b,timeVec]=inverse_read_stc(fn);
		end;
		cd(pdir);
		stc=mean(stc,3);

		base_idx=find(timeVec<0);

		z=etc_z(stc,base_idx);

		inverse_write_stc(z,v,a,b,sprintf('%s-%s.stc',output_stem{f_idx},hemi));
        end;
        %fn=sprintf('%s_%s-vol.stc',output_stem,erf_all(trig_idx).trig_str);
        %fprintf('\tsaving [%s]...\n',fn);
        %inverse_write_stc(X_dspm,[0:size(X_dspm,1)-1],min(erf_all(trig_idx).timeVec),t0,fn);

end;
cd(pdir);
fprintf('Done!\n');

