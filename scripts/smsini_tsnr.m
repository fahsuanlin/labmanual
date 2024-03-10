close all; clear all;

stc_file_stem={
    'smsini_mb_run_1_acc';
    'smsini_mb_run_2_acc';
    };

for file_idx=1:length(stc_file_stem)
    for hemi=1:2
        switch hemi
            case 1
                hemi_str='lh';
            case 2
                hemi_str='rh';
        end;

        fn=sprintf('%s-%s.stc',stc_file_stem{file_idx},hemi_str);
        [stc,v{hemi}]=inverse_read_stc(fn);

    
        [tsnr{file_idx,hemi}]=etc_tsnr('stc',stc,'exclude_time',[1:10]);

        fn_out=sprintf('%s_tsnr-%s.stc',stc_file_stem{file_idx},hemi_str);

        inverse_write_stc(repmat(tsnr{file_idx,hemi}(:),[1 5]),v{hemi},0,100,fn_out);
    end;

    for hemi=1:2
        switch hemi
            case 1
                hemi_str='lh';
            case 2
                hemi_str='rh';
        end;

        etc_render_fsbrain('subject','011624','overlay_vertex',v{hemi},'overlay_value',tsnr{file_idx,hemi}(:),'overlay_threshold',[50 500]);


        %generate image outputs
        for view_idx=1:2
            switch(hemi_str)
                case 'lh'
                    if(view_idx==1)
                        view(90,0);
                        view_str='med';
                    elseif(view_idx==2)
                        view(-90,0);
                        view_str='lat';
                    end;
                case 'rh'
                    if(view_idx==1)
                        view(-90,0);
                        view_str='med';
                    elseif(view_idx==2)
                        view(90,0);
                        view_str='lat';
                    end;
            end;
            hgexport(gcf,sprintf('%s_%s_%s.png',stc_file_stem{file_idx},hemi_str,view_str), hgexport('factorystyle'),'Format','png');
        end;
        close;
    end;
end;

save smsini_tsnr.mat tsnr;