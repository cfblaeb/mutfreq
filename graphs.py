__author__ = 'laeb'

import numpy as np
import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.widgets import Slider
from matplotlib import gridspec

import seaborn as sns


def heatplot(data, ax1, ax2):
	#data = pd.read_csv('combined.tsv', sep='\t', header=[0, 1], index_col=0, decimal=',')
	xstart = 7
	xend = 114
	# indicate everything we tested (0-1,3-4,6-7)
	# indicate our positive (1,4,7)
	# indicate their positive (3-5)
	# perhaps indicate wt (-1)
	our_test = (0, 1, 3, 4, 6, 7)
	our_positive = (1, 4, 7)
	their_positive = (3, 4, 5)

	ax2.spines['right'].set_visible(False)
	ax2.spines['top'].set_visible(False)
	ax2.yaxis.set_ticks_position('left')
	ax2.xaxis.set_ticks_position('bottom')

	ax2.set_yticks(np.arange(0.5, len(data.index), 1))
	ax2.set_yticklabels(data.index, ha='center', size=10)
	ax2.set_xticks(np.arange(0.5, len(data.columns), 1))
	ax2.set_xticklabels([x[1] for x in data.columns])

	our_plot_data = []
	their_plot_data = []

	for x, column in enumerate(data):
		our_plot_data.append(sum(list(data[column]).count(i) for i in our_positive) / sum(list(data[column]).count(i) for i in our_test))
		their_plot_data.append(sum(list(data[column]).count(i) for i in their_positive) / 20)
		for y, value in enumerate(data[column]):
			ax2.add_patch(mpl.patches.Rectangle((x, y), 1, 1, facecolor="None", edgecolor='darkgrey', alpha=0.2))
			if xstart < x < xend:
				if value in our_test and value not in our_positive:
					ax2.add_patch(mpl.patches.Rectangle((x, y), 1, 1, color="blue", alpha=0.6))
				if value == -1:
					ax2.add_patch(mpl.patches.Rectangle((x, y), 1, 1, facecolor="None", edgecolor='gold', alpha=0.6, zorder=999))
				if value in our_positive:
					#ax2.text(x + 0.5, y + 0.5, "O", size=12, family='monospace', ha='center', va='center')
					ax2.add_patch(mpl.patches.Rectangle((x, y), 1, 1, color="lightgreen"))
				if value in their_positive:
					ax2.text(x + 0.5, y + 0.5, "x", size=10, color='red', family='monospace', ha='center', va='center',)

			#ax2.text(x + 0.5,y + 0.5, value, size=8, family='monospace', ha='center', va='center', bbox=dict(facecolor='red', alpha=0.5))
			#ax2.add_patch(mpl.patches.Rectangle((x - .5, y - .5), 0.5, 0.5, facecolor="grey"))
	ax1.plot(np.arange(0.5, len(data.columns), 1), our_plot_data)
	ax1.plot(np.arange(0.5, len(data.columns), 1), their_plot_data)
	ax2.set_xlim(xstart + 1, xend)
	ax2.set_ylim(0, 20)


def enrichdata(control_freq, sorted_freq):
	for column in control_freq:  # remove WT
		control_freq[column][column[1]] = 0
		sorted_freq[column][column[1]] = 0

	enrichment = (sorted_freq / control_freq).loc[~sorted1.index.isin(['^', 'X', '+', '*'])].replace(np.inf, 0).replace(np.nan, 0)
	# lower outliers (above threshold in sorted, but virtually undetected in control)
	vals = enrichment.values.flatten()
	outlier_border = np.percentile(vals[vals > 0], 95)
	enrichment[enrichment > outlier_border] = outlier_border
	enrichment = enrichment / outlier_border
	return enrichment


def barchart(enricheddata):
	rasmol = {'A':'#C8C8C8','R':'#145AFF','N':'#00DCDC','D':'#E60A0A','C':'#E6E600','Q':'#00DCDC','E':'#E60A0A','G':'#EBEBEB','H':'#8282D2','I':'#0F820F','L':'#0F820F','K':'#145AFF','M':'#E6E600','F':'#3232AA','P':'#DC9682','S':'#FA9600','T':'#FA9600','W':'#B45AB4','Y':'#3232AA','V':'#0F820F'}
	rasmol_cmap = mpl.colors.ListedColormap([rasmol[x] for x in enricheddata.index])
	ax = enricheddata.T.plot(kind='bar', stacked=True, colormap=rasmol_cmap)
	for col in enricheddata:
		cum = 0
		for idx, val in enricheddata[col].iteritems():
			if val > 0:
				ax.annotate(idx, (col[0], cum + val / 2), ha='center', va='center')
			cum += val
	return ax

control1 = pd.read_csv('data/control_2_1.tsv', sep='\t', header=[0, 1], index_col=0, decimal=',')
sorted1 = pd.read_csv('data/sorted_2_1.tsv', sep='\t', header=[0, 1], index_col=0, decimal=',')

control2 = pd.read_csv('data/control_2_2.tsv', sep='\t', header=[0, 1], index_col=0, decimal=',')
sorted2 = pd.read_csv('data/sorted_2_2.tsv', sep='\t', header=[0, 1], index_col=0, decimal=',')

# The foldx data starts after index 7 and stops at 114
foldx = pd.read_csv('data/foldx.tsv', sep='\t', header=[0, 1], index_col=0, decimal=',')

miseq_noise_threshold = 0.003   # Miseq measured mutations present at frequencies less than this may be artifacts so we remove them
enrichment_threshold = 0.25		# Enrichment should ideally be exactly 1 for the majority of true positives. Anything less indicates less enrichment or more likely simply noise (e.g. something big in pre-sorting, will also be big post-sorting, it just wont be enriched)
foldx_threshold = 5

control_freq1 = control1 / control1.sum()
sorted_freq1 = sorted1 / sorted1.sum()
sorted_freq1[sorted_freq1 < miseq_noise_threshold] = 0  # filter out noise

control_freq2 = control2 / control2.sum()
sorted_freq2 = sorted2 / sorted2.sum()
sorted_freq2[sorted_freq2 < miseq_noise_threshold] = 0  # filter out noise

enrich1 = enrichdata(control_freq1, sorted_freq1)
enrich2 = enrichdata(control_freq2, sorted_freq2)
enrich1[enrich1 < enrichment_threshold] = 0  # filter out noise
enrich2[enrich2 < enrichment_threshold] = 0  # filter out noise

sns.set(style="white", context="talk")
#barchart(enrich1)
#barchart(enrich2)
#sns.heatmap(enrich1, cmap='YlOrRd')
#sns.heatmap(enrich2, cmap='Purples', alpha=0.5)
#plt.show()


def combinemap(control_freq1, control_freq2, sorted_freq1, sorted_freq2, enrich1, enrich2, foldx, foldx_threshold):
	#convert enrich1, 2 and foldx to
	# 0 = no test
	# 1 = neg
	# 2 = pos

	valuemap = enrich1.copy()
	valuemap[(control_freq1 > 0) | (control_freq2 > 0) | (sorted_freq1 > 0) | (sorted_freq2 > 0)] = 1
	valuemap[(enrich1 > 0) | (enrich2 > 0)] = 2
	#sns.heatmap(valuemap)

	foldxmap = foldx.copy()
	foldxmap.fillna(0, inplace=True)
	foldxmap[foldx < foldx_threshold] = 1
	foldxmap[foldx >= foldx_threshold] = 2
	foldxmap[foldx.isnull()] = 0
	#sns.heatmap(foldxmap)
	#plt.show()

	#then combine in -1 to 8 system and feed to heatplot
	"""
	-1=WT
		They predicted	Negative	Positive	We tested	Negative	Positive
	0	Yes				X						Yes			X
	1	Yes				X						Yes						X
	2	Yes				X						No
	3	Yes							X			Yes			X
	4	Yes							X			Yes						X
	5	Yes							X			No
	6	No										Yes			X
	7	No										Yes						X
	8	No										No
	"""
	combmap = enrich1.copy()
	combmap[(foldxmap == 1) & (valuemap == 1)] = 0
	combmap[(foldxmap == 1) & (valuemap == 2)] = 1
	combmap[(foldxmap == 1) & (valuemap == 0)] = 2
	combmap[(foldxmap == 2) & (valuemap == 1)] = 3
	combmap[(foldxmap == 2) & (valuemap == 2)] = 4
	combmap[(foldxmap == 2) & (valuemap == 0)] = 5
	combmap[(foldxmap == 0) & (valuemap == 1)] = 6
	combmap[(foldxmap == 0) & (valuemap == 2)] = 7
	combmap[(foldxmap == 0) & (valuemap == 0)] = 8

	#change wt to -1
	for column in combmap:  # remove WT
		combmap[column][column[1]] = -1
	return combmap

combinedmap = combinemap(control_freq1, control_freq2, sorted_freq1, sorted_freq2, enrich1, enrich2, foldx, foldx_threshold)

def fancyplot():
	styles = ['dark_background', 'grayscale', 'ggplot', 'fivethirtyeight', 'bmh']
	styles_choice = 0
	plt.style.use(styles[styles_choice])

	fig = plt.figure(figsize=(20, 8))

	gs = gridspec.GridSpec(2, 1, 0.05, 0.25, 0.95, 0.95, 0.05, 0.02)
	#fig, axes = plt.subplots(nrows=2, ncols=1, sharex=True, figsize=(20, 6))

	ax1 = fig.add_subplot(gs[0])
	ax2 = fig.add_subplot(gs[1], sharex=ax1)
	#plt.setp(ax2.get_xticklabels(), visible=False)
	heatplot(combinedmap, ax1, ax2)
	#sns.heatmap(combinedmap)
	#plt.show()
	#now make controls...foldx threshold, enrichment_threshold, perhaps miseq_threshold
	axcolor = 'lightgoldenrodyellow'
	axen_th = plt.axes([0.25, 0.1, 0.65, 0.03], axisbg=axcolor)
	axfx_th  = plt.axes([0.25, 0.15, 0.65, 0.03], axisbg=axcolor)

	en_th = Slider(axen_th, 'enrichment_threshold', 0.0, 1.0, valinit=enrichment_threshold)
	fx_th = Slider(axfx_th, 'foldx_threshold', -10, 40, valinit=foldx_threshold)

	def update(val):
		global old_amp, enrich1, enrich2
		amp = en_th.val
		if amp != old_amp:
			enrich1 = enrichdata(control_freq1, sorted_freq1)
			enrich2 = enrichdata(control_freq2, sorted_freq2)
			enrich1[enrich1 < amp] = 0  # filter out noise
			enrich2[enrich2 < amp] = 0  # filter out noise
			old_amp = amp
		freq = fx_th.val
		ax1.cla()
		ax2.cla()
		combinedmap = combinemap(control_freq1, control_freq2, sorted_freq1, sorted_freq2, enrich1, enrich2, foldx, freq)
		heatplot(combinedmap, ax1, ax2)
		fig.canvas.draw_idle()
	en_th.on_changed(update)
	fx_th.on_changed(update)

	plt.show()

old_amp = enrichment_threshold
#fancyplot()

#Det følgende er til at lave sens vs pres graf
xthres = np.arange(-5, 35, 1)
sensitivity = []  # how many percent of true positives (1,4) are found (4)
precision = []    # how many percent of positives (3,4) are true positives (4)

for xt in xthres:
	df = combinemap(control_freq1, control_freq2, sorted_freq1, sorted_freq2, enrich1, enrich2, foldx, xt)
	one = np.sum(~np.isnan(df[df == 1].values.flatten()))
	three = np.sum(~np.isnan(df[df == 3].values.flatten()))
	four = np.sum(~np.isnan(df[df == 4].values.flatten()))
	sensitivity.append(four/(four+one))
	precision.append(four/(four+three))

plt.plot(xthres, sensitivity)
plt.plot(xthres, precision)
plt.grid()
plt.show()
"""
	-1=WT
		They predicted	Negative	Positive	We tested	Negative	Positive
	0	Yes				X						Yes			X
	1	Yes				X						Yes						X
	2	Yes				X						No
	3	Yes							X			Yes			X
	4	Yes							X			Yes						X
	5	Yes							X			No
	6	No										Yes			X
	7	No										Yes						X
	8	No										No
"""
