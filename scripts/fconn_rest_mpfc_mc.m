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


mc_mat={
	'mc_037.mat';
        'mc_037.mat';
        'mc_037.mat';
        'mc_037.mat';
        'mc_037.mat';
        'mc_037.mat';
        'mc_037.mat';
        'mc_037.mat';
        'mc_037.mat';
        'mc_037.mat';
        'mc_037.mat';
};

fstem_confound={
	'regressor_wm_ventrical_037';
        'regressor_wm_ventrical_037';
        'regressor_wm_ventrical_037';
        'regressor_wm_ventrical_037';
        'regressor_wm_ventrical_037';
        'regressor_wm_ventrical_037';
        'regressor_wm_ventrical_037';
        'regressor_wm_ventrical_037';
        'regressor_wm_ventrical_037';
        'regressor_wm_ventrical_037';
        'regressor_wm_ventrical_037';
};


seed_str={
	'seed_mpfc_label';
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
                
                STC=cat(1,STC,stc{hemi_idx});
                flag_fe=1;
            else
                flag_fe=0;
            end;
        end;
        if(flag_fe)
            
            valid_subj_idx=cat(1,valid_subj_idx,d_idx);
            
            %fn=sprintf('%s/analysis/%s_regressors.mat',data_path,subject{d_idx});
	    %fn=sprintf('%s/analysis/%s_%s_regressors.mat',data_path,fstem_confound{f_idx},subject{d_idx});
            fn=sprintf('%s/%s/analysis/%s.mat',data_path,subject{d_idx},fstem_confound{f_idx});

            if(exist(fn))
                D_reg=[];
		fprintf('confound file [%s] found ...\n',fn);
                
		fconn{d_idx}.file_confound{f_idx}=fn;

		load(fn);
                %D_reg(:,1)=regressor_ventricle(1:end-1);
                %D_reg(:,2)=regressor_wm(1:end-1);

		regressor_ventricle=regressor_ventricle(1:size(STC,2));
		regressor_wm=regressor_wm(1:size(STC,2));
		D_reg(:,1)=regressor_ventricle;
		D_reg(:,2)=regressor_wm;

                D_reg(1:n_dummy,:)=[];
                D_reg(end-n_dummy+1:end,:)=[];
            else
                D_reg=[];
            end;

	    %motion regressors
            %fn=sprintf('%s/%s/Preprocessed/Movement_Regressors.txt',data_path,subject{d_idx});
	    fn=mc_mat{f_idx};
            if(exist(fn))
		fprintf('confound file [%s] found ...\n',fn);
		try
		%D_mc=textread(fn);
		D_mc=load(fn);
		D_mc=D_mc.mc;

		D_mc=D_mc(1:size(STC,2));
		%D_mc(end,:)=[];
		%D_mc(end,:)=[];
		catch
			fprintf('error in reading motion regressor!\n');
			D_mc=[];
		end;
	    else
		D_mc=[];
	    end;

	    try
	    	D_reg=cat(2,D_reg,D_mc);
            catch
		fprintf('error in concatenating WM/ventricle regressors [%s] and motion regressors [%s]!\n',mat2str(size(D_reg)),mat2str(size(D_mc)));
	    end;
            %remove global mean
            D=ones(size(STC,2),1);
	    D(:,2)=[0:size(D,1)-1]./(size(D,1)-1);
            if(~isempty(D_reg))
                D=cat(2,D,D_reg);
            end;
            if(flag_gavg);
                D=cat(2,D,mean(STC,1)');
            end;
            STC=(STC'-D*(inv(D'*D)*D'*STC')).';
            
            stc_hemi{1}=STC(1:length(v{1}),:);
            stc_hemi{2}=STC(length(v{1})+1:end,:);


	    hippo_reg=[];
	    for r_idx=1:1
	    	%file_regressor=sprintf('%s/%s/analysis/%s.mat',data_path,subject{d_idx},fstem_reg{f_idx,r_idx});

	    	%if(exist(file_regressor))
		%	fprintf('regressor file [%s] found...\n',file_regressor);

	        %        fconn{d_idx}.fconn{f_idx,r_idx}.file_regressor=file_regressor;

            	%	load(file_regressor);

                        %regressor_hippo=regressor_hippo(:)-(D*(inv(D'*D)*D'*regressor_hippo(:)));

		%	hippo_reg(:,r_idx)=regressor_hippo(1:end-1);

			%regressor_hippo=STC(9884,:)';

            		for hemi_idx=1:2
            			%%% do your analysis for each subject here ....
				%tmp=regressor_hippo(1:end-1);
				%tmp=regressor_hippo_full(:,1:end-1);
				%tmp=regressor_hippo;
				%tmp=tmp(:)-D*inv(D'*D)*D'*tmp(:);

				switch hemi_idx
				case 1
					ll=inverse_read_label('fconn_hcp_hippo_aseg_mPFC-lh.label');
					ll_offset=0;
				case 2
                                        ll=inverse_read_label('fconn_hcp_hippo_aseg_mPFC-rh.label');
					ll_offset=10242;
				end;
				[dummy,ll_now]=intersect([0:10241],ll);
				tmp=mean(STC(ll_offset+ll_now,:),1);
				
				%tmp=mean(tmp.'-D*inv(D'*D)*D'*tmp.',2);
	    			rr=etc_corrcoef(tmp(:),stc_hemi{hemi_idx}.');
            			fconn{d_idx}.fconn{f_idx,r_idx}.fconn(:,hemi_idx)=0.5.*log((1+rr)./(1-rr))./(1/sqrt(size(STC,2)/(2/TR*2.34)-3));   
				%%% end of subject-wise analysis.....
			end;		
	    	%end;
	    end;

	    if(~isempty(hippo_reg))
		%average of left and right hippocampus
	    	hippo_reg_avg=mean(hippo_reg,2);

		fconn{d_idx}.fconn{f_idx,r_idx+1}.fille_regressor='average';
		for hemi_idx=1:2
                                %%% do your analysis for each subject here ....
                                tmp=hippo_reg_avg(1:end);
				tmp=tmp(:)-D*inv(D'*D)*D'*tmp(:);
				rr=etc_corrcoef(tmp(:),stc_hemi{hemi_idx}.');
                                fconn{d_idx}.fconn{f_idx,r_idx+1}.fconn(:,hemi_idx)=0.5.*log((1+rr)./(1-rr))./(1/sqrt(size(STC,2)/(2/TR*2.34)-3));
                                %%% end of subject-wise analysis.....
            	end;
	    end;
        end;
    end;

    %aggregate fconn across cortical depths; One file for each seed ROI with stc indices across depths.   
    for r_idx=1:size(fconn{d_idx}.fconn,2)	
    	for hemi_idx=1:2
		switch hemi_idx
		case 1
			hemi='lh';
		case 2
			hemi='rh';
		end;

		fn=sprintf('fconn_%s-%s.stc',seed_str{r_idx},hemi);
		tmp=[];
		for f_idx=1:size(fconn{d_idx}.fconn,1)
			tmp(:,f_idx)=fconn{d_idx}.fconn{f_idx,r_idx}.fconn(:,hemi_idx);
		end;
		inverse_write_stc(tmp,[0:10241],0,100,fn);

	end;
    end;
 
    fprintf('\n');

end;
save fconn_rest_mpfc_label_mc.mat fconn
