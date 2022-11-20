#!/usr/bin/env python

import numpy as np
import matplotlib.pyplot as plt
import nibabel as nib

import rpy2.robjects as ro
from rpy2.robjects import numpy2ri
from rpy2.robjects.packages import importr

## suppress console because of weird permission around r
#import warnings
#from rpy2.rinterface import RRuntimeWarning
#warnings.filterwarnings("ignore", category=RRuntimeWarning)
#import warnings
#warnings.filterwarnings('ignore')


# load r functions
stats = importr('stats')
base = importr('base')
utils = importr('utils')

ggplot2 = importr('ggplot2')
dplyr = importr('dplyr')
scales = importr('scales')

ggseg = importr('ggseg')
ggsegGlasser = importr('ggsegGlasser')



scriptwd = os.path.dirname(os.path.dirname(os.path.realpath(__file__)))
be_template = "%s/ants_tbss/ants_oasis_template_ras/T_template0.nii.gz" % scriptwd

numpy2ri.activate()


def concert_aseg(atlasdata, dataname, labels = None):
	if labels is None:
		labels = np.genfromtxt("aseg_labels.txt", dtype = str)
		regions = np.array(labels)
		hemi = np.repeat('midline', len(labels))
		for l, label in enumerate(labels):
			if label.startswith("Left-"):
				hemi[l] = 'left'
				temp = label.replace("Left-","")
				temp = temp.lower()
				temp = temp.replace("-"," ")
				regions[l] = temp
			if label.startswith("Right-"):
				hemi[l] = 'right'
				temp = label.replace("Right-","")
				temp = temp.lower()
				temp = temp.replace("-"," ")
				regions[l] = temp
	ratlasdata = ro.r.tibble(hemi = hemi, region = regions, label = labels, value = atlasdata)
	ratlasdata.names = ['hemi', 'region', 'label', dataname]
	return(ratlasdata)

def convert_glasser(atlasdata, dataname,
			lannot = '/media/tris/PortableSSD/RESULTS_WRITEUP_SCCA/RESULTS/SURFACE/HCP/lh.hcp-mmp-b-fix.annot',
			rannot = '/media/tris/PortableSSD/RESULTS_WRITEUP_SCCA/RESULTS/SURFACE/HCP/rh.hcp-mmp-b-fix.annot'):
	labels, ctab, names = nib.freesurfer.read_annot(lannot)
	lhnames = np.array(np.array(names)[np.unique(labels)[1:]], str)
	lhnames = np.array([("lh_" + name.replace("_ROI", "")) for name in lhnames])
	lhhemi = np.repeat("left", len(lhnames))
	labels, ctab, names = nib.freesurfer.read_annot(rannot)
	rhnames = np.array(np.array(names)[np.unique(labels)[1:]], str)
	rhnames = np.array([("rh_" + name.replace("_ROI", "")) for name in rhnames])
	rhhemi = np.repeat("right", len(rhnames))
	regions = np.concatenate((lhnames, rhnames))
	hemi = np.concatenate((lhhemi, rhhemi))
	labels = np.array(regions)
	regions = np.array([name[5:] for name in regions]) 
	ratlasdata = ro.r.tibble(hemi = hemi, region = regions, label = labels, value = atlasdata)
	ratlasdata.names = ['hemi', 'region', 'label', dataname]
	return(ratlasdata)








#colours = c("red","yellow","green","lightblue","darkblue","lightblue","green","yellow","red"),
#                         values = c(1.0,0.8,0.6,0.4,0.2,0,0.2,0.4,0.6,0.8,1.0)

out_components = np.load('outcomponents.npy')
metric = 'loading'
for i in range(len(out_components)):
	ro.globalenv['ratlasdata'] = convert_glasser(atlasdata = out_components[i,:358], dataname = metric)

	fig=ro.r("""
	  ggplot(ratlasdata) +
	  geom_brain(atlas = glasser, 
		          position = position_brain(hemi ~ side),
		          aes(fill = %s)) +
		          scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0, guide = "colourbar") +
		          theme_void() +
		          labs(title = "Cortical Comp [%d]")
	""" % (metric, int(i+1)))
	ro.r.ggsave(
	  "cort_output_%d.png" % int(i+1),
	  fig,
	  width = 4,
	  height = 4,
	  dpi = 1200
	)

	fig=ro.r("""
	  ggplot(ratlasdata) +
	  geom_brain(atlas = glasser, 
		          position = position_brain(hemi ~ side),
		          aes(fill = %s)) +
		          scale_fill_viridis_c(option = "cividis", direction = -1) +
		          theme_void() +
		          labs(title = "Cortical Comp [%d]")
	""" % (metric, int(i+1)))
	ro.r.ggsave(
	  "cort_output_%d.png" % int(i+1),
	  fig,
	  width = 4,
	  height = 4,
	  dpi = 1200
	)



	ro.globalenv['ratlasdata'] = concert_aseg(atlasdata = out_components[i,358:], dataname = metric)
	fig=ro.r("""
	  ggplot(ratlasdata) +
	  geom_brain(atlas = aseg, 
		          aes(fill = %s)) +
		          scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0, guide = "colourbar") +
		          theme_void() +
		          labs(title = "Subcortical Comp [%d]")
	""" % (metric, int(i+1)))
	ro.r.ggsave(
	  "subcort_output_%d.png" % int(i+1),
	  fig,
	  width = 4,
	  height = 4,
	  dpi = 1200
	)

#limits=c(-0.31, 0.31), oob=squish
p=ro.r("""
ggplot(ratlasdata) +
  geom_brain(atlas = glasser, 
             position = position_brain(hemi ~ side),
             aes(fill = p)) +
  scale_fill_viridis_c(option = "cividis", direction = -1) +
  theme_void() +
  labs(title = "My awesome title", 
       subtitle = "of a brain atlas plot",
       caption = "I'm pretty happy about this!")
""")
