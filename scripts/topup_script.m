close all; clear all;

dicom_path_pa={
	'012';
};

dicom_path_ap={
	'009';
	'010';
	'011';
	'020';
	'021';
	'022';
	'023';
};

dicom_path_root='/space_lin2/fhlin/7t_music_skku/LAM_AUD_SHL/unpack/bold';
dicom_file_name='fmcprstc.nii.gz';

[dummy,dicom_file_name_stem]=fileparts(dicom_file_name);


for f_idx=1:length(dicom_path_pa)
	cmd=sprintf('!fslroi %s/%s/%s %s/%s/PA_b0 0 1',dicom_path_root, dicom_path_pa{f_idx}, dicom_file_name, dicom_path_root, dicom_path_pa{f_idx});
	eval(cmd);
end;

for f_idx=1:length(dicom_path_ap)
        cmd=sprintf('!fslroi %s/%s/%s %s/%s/AP_b0 0 1',dicom_path_root, dicom_path_ap{f_idx}, dicom_file_name, dicom_path_root, dicom_path_ap{f_idx});
	eval(cmd);

	cmd=sprintf('!fslmerge -t %s/%s/AP_PA_%s %s/%s/AP_b0.nii.gz %s/%s/PA_b0.nii.gz',dicom_path_root, dicom_path_ap{f_idx},dicom_file_name, dicom_path_root, dicom_path_ap{f_idx}, dicom_path_root, dicom_path_pa{1});
	eval(cmd);

	cmd=sprintf('!topup --imain=%s/%s/AP_PA_%s --datain=acqp.txt --subsamp=1 --config=b02b0.cnf --out=%s/%s/AP_PA_fmcprstc_topup --fout=%s/%s/AP_PA_fmcprstc_topup_field --iout=%s/%s/AP_PA_fmcprst_upwarp_b0 -v',dicom_path_root, dicom_path_ap{f_idx},dicom_file_name,dicom_path_root, dicom_path_ap{f_idx}, dicom_path_root, dicom_path_ap{f_idx},dicom_path_root, dicom_path_ap{f_idx});

	eval(cmd);
	cmd=sprintf('!applytopup --imain=%s/%s/%s  --inindex=1 -a acqp.txt --topup=%s/%s/AP_PA_fmcprstc_topup --out=%s/%s/fmcprstc_topup.nii.gz --method=jac', dicom_path_root, dicom_path_ap{f_idx},dicom_file_name, dicom_path_root, dicom_path_ap{f_idx}, dicom_path_root, dicom_path_ap{f_idx});

	eval(cmd);
        fprintf('[%s] OK!\n',dicom_path_ap{f_idx});

end;


