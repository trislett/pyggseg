#!/usr/bin/env python

import os

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


# get the working directory
scriptwd = os.path.dirname(os.path.dirname(os.path.realpath(__file__)))

def convert_aseg(atlasdata, dataname, label_file = "%s/static/aseg_labels.txt" % scriptwd):
	"""
	Converts the freesurfer automatic segmentation (aseg) to an R object that can be read by ggseg.
	Unfortunately, the accumbens is out of view. This isn't necessary IMO. Maybe I'll fix it for them later...
	
	Parameters
	----------
	atlasdata : array
		The data files that corresponds to the atlas labels
	dataname : str
		The output name of variables. e.g., p-value, rho, etc.
	label_file : str
		The location of label file to import. The length of the label_file must match the atlasdata (and in the same order).
	Returns
	-------
	ratlasdata : robj
		A tibble r object that is be passed to R global environment as 'ratlasdata'.
	"""
	numpy2ri.activate()
	labels = np.genfromtxt(label_file, dtype = str)
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
	ratlasdata.dtype.names = ('hemi', 'region', 'label', dataname)
	ro.globalenv['ratlasdata'] = ratlasdata
	numpy2ri.deactivate()
	return(ratlasdata)

def convert_glasser(atlasdata, dataname,
			lannot = "%s/static/lh.hcp-mmp-b-fix.annot" % scriptwd,
			rannot = "%s/static/rh.hcp-mmp-b-fix.annot" % scriptwd):
	"""
	Converts a modified version of the HCP_MMP1 atlas to an R object that can be read by ggseg.
	
	Parameters
	----------
	atlasdata : array
		The data files that corresponds to the atlas labels
	dataname : str
		The output name of variables. e.g., p-value, rho, etc.
	label_file : str
		The location of label file to import. The length of the label_file must match the atlasdata (and in the same order).
	Returns
	-------
	ratlasdata : robj
		A tibble r object that is be passed to R global environment as 'ratlasdata'.
	"""
	numpy2ri.activate()
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
	ratlasdata.dtype.names = ('hemi', 'region', 'label', dataname)
	ro.globalenv['ratlasdata'] = ratlasdata
	numpy2ri.deactivate()
	return(ratlasdata)

def run_test_hcp(output_file = "glasser_test_output.png", maxscale = None, title = None):
	values = np.random.random(358)
	values[179:] = values[179:] * -1
	metric = "pvalues"
	if maxscale is None:
		maxscale = np.ceil(np.abs(values).max() * 1000.) / 1000.
		print("%1.3f" % maxscale)
	convert_glasser(atlasdata = values, dataname = metric)
	cmd = ("""
	ggplot(ratlasdata) +
	geom_brain(atlas = glasser, 
					position = position_brain(hemi ~ side),
					aes(fill = %s)) +
					scale_fill_gradientn(colors = c("#00008C", "#2234A8", "#4467C4", "#659BDF", "#87CEFB", "white", "#ffec19", "#ffc100", "#ff9800", "#ff5607", "#f6412d"), limits = c(%1.3f, %1.3f), guide = "colourbar") +
					theme_void()""" % (metric, -maxscale, maxscale))
	if title is not None:
		cmd = cmd + """ +\nlabs(title = "%s")""" % title
	fig = ro.r(cmd)
	ro.r.ggsave(output_file, fig, width = 6, height = 4, dpi = 1200)

def run_test_aseg(output_file = "aseg_test_output.png", maxscale = None, title = None):
	values = np.random.random(14)
	values[7:] = values[7:] * -1
	metric = "pvalues"
	if maxscale is None:
		maxscale = np.ceil(np.abs(values).max() * 1000.) / 1000.
		print("%1.3f" % maxscale)
	convert_aseg(atlasdata = values, dataname = metric)
	cmd = ("""
	ggplot(ratlasdata) +
	geom_brain(atlas = aseg, 
					side = "coronal",
					aes(fill = %s)) +
					scale_fill_gradientn(colors = c("#00008C", "#2234A8", "#4467C4", "#659BDF", "#87CEFB", "white", "#ffec19", "#ffc100", "#ff9800", "#ff5607", "#f6412d"), limits = c(%1.3f, %1.3f), guide = "colourbar") +
					theme_void()""" % (metric, -maxscale, maxscale))
	if title is not None:
		cmd = cmd + """ +\nlabs(title = "%s")""" % title
	fig = ro.r(cmd)
	ro.r.ggsave(output_file, fig, width = 4, height = 3, dpi = 1200)



