close all; clear all;

file_stem={
	'tfMRI_WM_surf_soa_glm_h01_effect';
        'tfMRI_WM_surf_soa_glm_h02_effect';
        'tfMRI_WM_surf_soa_glm_h03_effect';
};

data_path='/space_lin1/hcp';

output_stem={
	'tfMRI_WM_soa_glm_h01';
        'tfMRI_WM_soa_glm_h02';
        'tfMRI_WM_soa_glm_h03';
};


subject='';
d=textread('subject_list_all.txt');
for d_idx=1:length(d)
    subject{d_idx}=num2str(d(d_idx));
end;
%subject=[];
%subject{1}='119732';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for f_idx=1:length(file_stem)
	for hemi=1:2
		switch hemi
		case 1
			hemi_str='lh';
		case 2
			hemi_str='rh';
		end;

		counter=1;
		for d_idx=1:length(subject)
    			fprintf('{%s}::[%s]...(%04d|%04d)....\r',file_stem{f_idx},subject{d_idx},d_idx,length(subject));

    			fn=sprintf('%s/%s/analysis/%s-%s.stc',data_path,subject{d_idx},file_stem{f_idx},hemi_str);
			if(exist(fn))
				[stc,v]=inverse_read_stc(fn);

				stc_all(f_idx).stc(:,counter,hemi)=stc(:,1);
				if(counter==1)
					tmp=stc(:,1);
					tmp2=stc(:,1).^2;
				else
					tmp=tmp+stc(:,1);
					tmp2=tmp2+stc(:,1).^2;
				end;
				counter=counter+1;
			end;
		end;
		fprintf('\n');
		stc_avg=tmp./counter;
		tmp2=tmp2./counter;
		stc_std=sqrt(tmp2-stc_avg.^2);
		stc_z=stc_avg./stc_std.*sqrt(counter-1);

		fn=sprintf('%s-z-%s.stc',output_stem{f_idx},hemi_str);
		inverse_write_stc(repmat(stc_z(:),[1 5]),v,0,100,fn);
	end;
end;

%save wm.mat stc_all
