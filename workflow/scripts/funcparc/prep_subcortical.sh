#!/bin/bash
# The subcortical fMRI timeseries need to get slightly smoothed (1.6mm FWHM) to match the 1.6mm average vertex spacing of the 59k_fs_LR mesh.

invol=$1
rois=$2
tempdir=$3
Sigma=$4
outvol=$5

NameOffMRI=${invol%.nii.gz*}
NameOffMRI=`basename $NameOffMRI`

mkdir -p $tempdir

echo "Processing subcortical fMRI timeseries to match surface data"
echo ""

{ echo "- Creating subject-roi subcortical cifti at same resolution as output"
echo -e "\twb_command -cifti-create-dense-timeseries $tempdir/${NameOffMRI}_temp_atlas.dtseries.nii -volume $invol $rois"
} | fmt -w 150

	wb_command -cifti-create-dense-timeseries $tempdir/${NameOffMRI}_temp_atlas.dtseries.nii -volume $invol $rois

{ echo "- Dilating out zeros"
echo -e "\twb_command -cifti-dilate $tempdir/${NameOffMRI}_temp_atlas.dtseries.nii COLUMN 0 30 $tempdir/${NameOffMRI}_temp_atlas_dilate.dtseries.nii"
} | fmt -w 150

	wb_command -cifti-dilate $tempdir/${NameOffMRI}_temp_atlas.dtseries.nii COLUMN 0 30 $tempdir/${NameOffMRI}_temp_atlas_dilate.dtseries.nii
	rm -f $tempdir/${NameOffMRI}_temp_atlas.dtseries.nii

{ echo "- Smoothing and resampling"
echo -e "\twb_command -cifti-smoothing $tempdir/${NameOffMRI}_temp_atlas_dilate.dtseries.nii 0 ${Sigma} COLUMN $tempdir/${NameOffMRI}_temp_atlas_smooth.dtseries.nii -fix-zeros-volume"
} | fmt -w 150

	wb_command -cifti-smoothing $tempdir/${NameOffMRI}_temp_atlas_dilate.dtseries.nii 0 ${Sigma} COLUMN $tempdir/${NameOffMRI}_temp_atlas_smooth.dtseries.nii -fix-zeros-volume

{ echo -e "\twb_command -cifti-dilate $tempdir/${NameOffMRI}_temp_atlas_smooth.dtseries.nii COLUMN 0 30 $tempdir/${NameOffMRI}_temp_atlas_dilate.dtseries.nii"
} | fmt -w 150

	wb_command -cifti-dilate $tempdir/${NameOffMRI}_temp_atlas_smooth.dtseries.nii COLUMN 0 30 $tempdir/${NameOffMRI}_temp_atlas_dilate.dtseries.nii
	rm -f $tempdir/${NameOffMRI}_temp_atlas_smooth.dtseries.nii

{ echo -e "\twb_command -cifti-separate $tempdir/${NameOffMRI}_temp_atlas_dilate.dtseries.nii COLUMN -volume-all $outvol"
} | fmt -w 150

	wb_command -cifti-separate $tempdir/${NameOffMRI}_temp_atlas_dilate.dtseries.nii COLUMN -volume-all $outvol

#	rm -rf $tempdir
