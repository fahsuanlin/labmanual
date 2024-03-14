close all; clear all;

subject={
    's025';
    's026';
    's031';
    's032';
    's034';
%    's041';
%    's045';
%    's047';
    };
    


cond_stem={
    'fmri_surf_soa_glm_h01_beta';
    'fmri_surf_soa_glm_h02_beta';
    'fmri_surf_soa_glm_h03_beta';
};

cond_output_stem={
    'h01';
    'h02';
    'h03';
};


cond_stem_str={
    'A';
    'V';
    'AV';
    };

root_dir='/Users/fhlin/workspace/seeg';

output_stem='average_surf_fdr_031324';

for cond_idx=1:length(cond_stem)
    z_both=[];
    for hemi_idx=1:2
        switch hemi_idx
            case 1
                hemi_str='lh';
            case 2
                hemi_str='rh';
        end;
        
        for subj_idx=1:length(subject)
            [dummy,v,a,b,timeVec]=inverse_read_stc(sprintf('%s/%s/fmri_analysis/%s-%s.stc',root_dir,subject{subj_idx},cond_stem{cond_idx},hemi_str));

            stc(:,subj_idx)=dummy(:,1);
        end;
        
        z=mean(stc,2)./std(stc,0,2).*sqrt(size(stc,2));
        
        inverse_write_stc(repmat(z(:),[1:5]),v,a,b,sprintf('%s_%s-%s.stc',output_stem,cond_output_stem{cond_idx},hemi_str));

        z_both=cat(1,z_both(:),z(:));
       
    end;
    p_both=1-normcdf(abs(z_both));

    [p_fdr5, p_masked] = fdr( p_both, 0.05);
    [p_fdr1, p_masked] = fdr( p_both, 0.01);
    
    %
    th_5p=norminv(1-p_fdr5);
    th_1p=norminv(1-p_fdr1);
    
    etc_render_fsbrain_stc({sprintf('%s_%s',output_stem,cond_output_stem{cond_idx})},[th_5p th_1p],'flag_overlay_pos_only',1);

    print('-dpng',sprintf('%s/average/fmri/%s_%s_average.png',root_dir,output_stem,cond_stem_str{cond_idx}));
    
    close

end;
