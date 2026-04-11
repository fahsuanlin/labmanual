close all; clear all;

data_path='/space_lin2/fhlin/sri_7t';

fstem={
	'2_fsaverage_fmcprstc_topup_037_pial';
        '2_fsaverage_fmcprstc_topup_037_gray09';
	'2_fsaverage_fmcprstc_topup_037_gray08';
        '2_fsaverage_fmcprstc_topup_037_gray07';
        '2_fsaverage_fmcprstc_topup_037_gray06';
        '2_fsaverage_fmcprstc_topup_037_gray05';
        '2_fsaverage_fmcprstc_topup_037_gray04';
        '2_fsaverage_fmcprstc_topup_037_gray03';
        '2_fsaverage_fmcprstc_topup_037_gray02';
        '2_fsaverage_fmcprstc_topup_037_gray01';
        '2_fsaverage_fmcprstc_topup_037_white';
    };



TR=2.5; %second

n_dummy=0;
flag_gavg=0;

subject{1}='hjlee';
%subject=[];
%subject{1}='119732';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for d_idx=1:length(subject)
    fprintf('[%s]...(%04d|%04d)....\r',subject{d_idx},d_idx,length(subject));
      
    fconn{d_idx}.subject=subject{d_idx};
 
    for f_idx=1:length(fstem)
	valid_subj_idx=[];

        roi=[];
        STC=[];

	ROI_avg=[];
	ROI_std=[];
        for hemi_idx=1:2
            switch hemi_idx
                case 1
                    hemi_str='lh';
                case 2
                    hemi_str='rh';
            end;
            
            fn=sprintf('%s/%s/analysis/%s_%s-%s.stc',data_path,subject{d_idx},subject{d_idx},fstem{f_idx},hemi_str);
            if(exist(fn))
		fprintf('data file [%s] found ...\n',fn);

		fconn{d_idx}.data_stem{f_idx}=fstem{f_idx};
                [stc{hemi_idx},v{hemi_idx},d0,d1,timeVec]=inverse_read_stc(fn);
                
                %remove dummy scans
                stc{hemi_idx}(:,1:n_dummy)=[];
                stc{hemi_idx}(:,end-n_dummy+1:end)=[];
    
		tsnr{hemi_idx}(:,f_idx)=etc_tsnr('stc',stc{hemi_idx});            
                flag_fe=1;
            else
                flag_fe=0;
            end;
        end;
    end;

    	for hemi_idx=1:2
		switch hemi_idx
		case 1
			hemi='lh';
		case 2
			hemi='rh';
		end;

		fn=sprintf('tsnr-%s.stc',hemi);
		inverse_write_stc(tsnr{hemi_idx},v{hemi_idx},0,100,fn);

	end;
 
    fprintf('\n');

end;
