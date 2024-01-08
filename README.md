![](https://github.com/fahsuanlin/labmanual/blob/master/images/background.png)

Welcome to the wiki of my lab! Here are lots of scripts for magnetic resonance imaging (MRI), magnetoencephalography (MEG), and electroencephalography (EEG) research. 

Goto the [wiki page](https://github.com/fahsuanlin/labmanual/wiki) for further details.

The web sites of our lab are [here](http://linbrainlab.org) and [here at Sunnybrook](https://sunnybrook.ca/research/content/?page=sri-groups-mbim-about).  

[Contact me](mailto:fhlin@sri.utoronto.ca) if there is any question/error.


## Structural MRI

- [FreeSurfer reconstruction](https://github.com/fahsuanlin/labmanual/wiki/0A.-FreeSurfer-reconstruction) for brain 3D models and anatomical segementation, labels, parameters .

- [FreeSurfer quick sheet](https://github.com/fahsuanlin/labmanual/wiki/36.-Quick-sheet-about-Freesurfer-reconstructions)

- [Diffusion tensor imaging (DTI) analysis](https://github.com/fahsuanlin/labmanual/wiki/17.-DTI-analysis)

## functional MRI (fMRI)
- [Stimuli presentation](https://github.com/fahsuanlin/labmanual/wiki/16.-Stimuli-presentation-and-response-collection-in-fMRI-by-Psychtoolbox) using Psychtoolbox
- [Unpack dicom files for analysis](https://github.com/fahsuanlin/labmanual/wiki/09.-unpack-dicom).
  
- [fMRI preprocessing](https://github.com/fahsuanlin/labmanual/wiki/10.-fMRI-data-preprocessing) includes motion correction, slice timing correction, and spatial smoothing.

- [Tasked fMRI analysis](https://github.com/fahsuanlin/labmanual/wiki/11.-fMRI-analysis) uses General Linear Modeling to identify brain regions correlated with a model of stimulatin, task performance, or other explicit (time-domain) models.

- [Resting-state fMRI analysis](https://github.com/fahsuanlin/labmanual/wiki/46.-Resting%E2%80%90state-fMRI-analysis) using 'seed-based' correlation to identify a network with correlated fMRI signal dynamics with a chosen brain location.

- [Inter-subject correlation analysis](https://github.com/fahsuanlin/labmanual/wiki/12.-fMRI-inter-subject-correlation-analysis) of fMRI using complex naturalistic stimuli.

- [SMS-InI analysis](https://github.com/fahsuanlin/labmanual/wiki/04.-SMS-InI-analysis)

- [Co-register between SMS-InI reconstructions and anatomical MRI](https://github.com/fahsuanlin/labmanual/wiki/38.-Registration-between-SMS%E2%80%90InI-and-anatomical-MRI)

- [Cortical-depth dependent fMRI analysis](https://github.com/fahsuanlin/labmanual/wiki/05.-layer-fMRI-analysis)


## EEG-MRI

- [Measurement setup](https://github.com/fahsuanlin/labmanual/wiki/20.-EEG-setup)

- [Data collection cheat sheet](https://github.com/fahsuanlin/labmanual/wiki/33.-EEG-fMRI-acquisition:-MRI-control)
  
- [EEG artifact suppression](https://github.com/fahsuanlin/labmanual/wiki/18.-Suppression-of-artifacts-in-EEG-collected-inside-MRI)

- [EEG evoked potentials and source analysis](https://github.com/fahsuanlin/labmanual/wiki/02.-EEG-analysis-stream)

## MEG 

- [MEG analysis](https://github.com/fahsuanlin/labmanual/wiki/03.-MEG-analysis-stream) for the data collected at Academia Sinica

## Stereotatic EEG (SEEG)

Some pages are from [fhlin_toolbox repo](https://github.com/fahsuanlin/fhlin_toolbox).

- [Register the location of SEEG electrodes using post-implantation MRI; Step I](https://github.com/fahsuanlin/fhlin_toolbox/wiki/SEEG:-register-electrodes-to-MRI)

- [Register the location of SEEG electrodes using post-implantation MRI; Step II](https://github.com/fahsuanlin/fhlin_toolbox/wiki/SEEG:-register-electrodes-to-MRI-(II))

- [Register the location of SEEG electrodes using post-implantation CT](https://github.com/fahsuanlin/fhlin_toolbox/wiki/SEEG:-register-electrodes-to-MRI-with-the-guidance-by-CT)

- [Visualize the registered electrodes](https://github.com/fahsuanlin/fhlin_toolbox/wiki/SEEG:-view-registered-electrodes)

- [Data analysis over electrodes](https://github.com/fahsuanlin/labmanual/wiki/07.-SEEG-electrode-analysis)

- [SEEG source analysis: individual](https://github.com/fahsuanlin/labmanual/wiki/06.-SEEG-source-modeling)

- [SEEG source analysis: group](https://github.com/fahsuanlin/labmanual/wiki/08.-SEEG-group-analysis)

## Transcranial magnetic stimulation

- [Image preparation for BrainSight](https://github.com/fahsuanlin/labmanual/wiki/23.-Image-processing-in-BrainSight)

- [Electric field modeling](https://github.com/fahsuanlin/labmanual/wiki/29.-TMS-E-field-modeling)

## Bruker 7T animal MRI

- [Data collection](https://github.com/fahsuanlin/labmanual/wiki/30.-Bruker-7T-scanning)

## Visualization

These are pages from [fhlin_toolbox repo](https://github.com/fahsuanlin/fhlin_toolbox).

- [Show 3D brain model for an individual](https://github.com/fahsuanlin/fhlin_toolbox/wiki/Render-brain)

- [Show 3D brain model and orthogonal slicse for an indvidual](https://github.com/fahsuanlin/fhlin_toolbox/wiki/Render-brain-and-check-brain-coordinates). This includes showing cortical labels created automatically by FreeSurfer (.annot file).

- [Show and create a region-of-interest on a 3D brain model](https://github.com/fahsuanlin/fhlin_toolbox/wiki/Render-brain-and-create-a-region-of-interest).

- [Show 3D brain model with cortical labels](https://github.com/fahsuanlin/fhlin_toolbox/wiki/Render-brain-labels-and-annotation). These labels are either from individual label files (*.label) or a collection of labels in an annotation file (.annot).

- [Show 3D brain model with overlay](https://github.com/fahsuanlin/fhlin_toolbox/wiki/Render-brain-with-an-overlay). The overlay is typically fMRI results, which is a static image here.

- [Show 3D brain model with dynamic overlay](https://github.com/fahsuanlin/fhlin_toolbox/wiki/Render-brain-with-a-dynamic-overlay). The overlay is typically fMRI results, which can be a static image or a collection of images over time.

- [Show or examine the spatial co-registration between structural MRI and fMRI](https://github.com/fahsuanlin/fhlin_toolbox/wiki/Render-or-register-brain-with-a-volume-overlay)


## Sample data

- [EEG-MRI: steady-state visual evoked potentials](https://github.com/fahsuanlin/labmanual/wiki/21.-Sample-data:-Steady-state-visual-potential)

- [fMRI: Physiological noise in healthy elderly and Alzheimer's disease patients](https://github.com/fahsuanlin/labmanual/wiki/22.-Sample-data:-physiological-noise-in-Alzheimer's-disease-patients-and-health-controls)

- [fMRI: perspective-taking effects probed by movies](https://github.com/fahsuanlin/labmanual/wiki/25.-Sample-data:-perspective-taking-fMRI-data)

- [EEG: Social perception in young adults and ADHD patients probed by movies](https://github.com/fahsuanlin/labmanual/wiki/31:-Sample-data:-social-interaction-EEG-data)

- [SEEG: music listening](https://github.com/fahsuanlin/labmanual/wiki/32:-Sample-data:-SEEG-recording-during-music-listening)

- [EEG-MRI: music listening](https://github.com/fahsuanlin/labmanual/wiki/39.-Sample-data:-EEG%E2%80%90MRI-during-music-listening)

- [MRI: Human Connectome Project young adults](https://github.com/fahsuanlin/labmanual/wiki/40:-Sample-data:-Human-Connectome-Project-(HCP)-fMRI)
  

## MRI coils
- [Coil housing printing](https://github.com/fahsuanlin/labmanual/wiki/43:-Coil-housing-printing)

- [Coil test rig for selective detuning](https://github.com/fahsuanlin/labmanual/wiki/24.-Coil:-TIM-4G-coil-test-rig)

- [3T receiver coil circuits](https://github.com/fahsuanlin/labmanual/wiki/41.-Coil:-3T-coil-schematics)

- [Impedance matching circuits for 3T and 7T](https://github.com/fahsuanlin/labmanual/wiki/43:-Coil-housing-printing)
  
- [3T receiver coil bench testing](https://github.com/fahsuanlin/labmanual/wiki/42:-Coil:-3T-coil-bench-testing)
