import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import sys
import scipy.cluster.vq as vq
import scipy.io as io
import sklearn.cluster as cluster
import matplotlib.cm as cm
import matplotlib.patheffects as pe
import argparse
from operator import itemgetter
from numpy import genfromtxt,zeros,linspace
from matplotlib.colors import BoundaryNorm
from matplotlib.colorbar import ColorbarBase,Colorbar
from scipy.signal import argrelextrema

figsize_spectrum = (25,8)

parser = argparse.ArgumentParser(description='Cluster eigenvectors and generate a pretty report.')
parser.add_argument('-t', metavar='TITLE', type=str, nargs='+', help='project title')
parser.add_argument('-x', metavar='MAT', type=str, nargs='+', help='Matlab matrix')
parser.add_argument('-y', metavar='TYPE', type=int, nargs='+', help='type of matlab input')
parser.add_argument('-e', metavar='VALUES', type=str, nargs='+', help='file of eigenvalues')
parser.add_argument('-v', metavar='VECTORS', type=str, nargs='+', help='file of eigenvectors')
parser.add_argument('-c', metavar='COORDS', type=str, nargs='+', help='file of coordinates')
parser.add_argument('-N', metavar='NUM_CLUSTERS', type=int, nargs='+', help='number of clusters')
parser.add_argument('-M', metavar='NUM_EIGENS', type=int, nargs='+', help='number of eigenvalues to consider')

args = parser.parse_args()

num_eigens = args.M[0]
num_centers = args.N[0]

matlab_type = 1
if args.y != None:
	matlab_type = args.y[0]

# my_palette_array = np.array([(256,256,256), (115,211,72), (200,79,140), (198,156,152), (143,211,150), (208,88,59), (145,147,203), (88,118,55), (80,49,92), (104,55,41), (104,79,184)]) / float(256)
my_palette_array = np.array([(1, 1, 1), (1, 0.5, 0), (0.25, 1, 1), (0.5, 0, 0), (1, 0, 0), (1, 1, 0), (0, 0, 1), (0.5, 0.5, 0.5), (0, 1, 0), (0, 0, 0.4), (0.75, 0.25, 0)])

def get_colormap(num_colors):
    # map cluster label to color (0, 1, 2, ...) -> (orange, blue, green, ...)  
    from matplotlib.colors import ListedColormap
    return ListedColormap(my_palette_array[0:(num_colors+1)]) 

# num_eigens = 10
# num_centers = 10

### read matrix
mat = io.loadmat(args.x[0])
spec_varname = 'SP'

### average spectrum plot
average_spectrum = np.mean(mat[spec_varname], axis=1)
len_spectrum = len(average_spectrum)
spec_x = range(1, len_spectrum+1)
average_spectrum_pairs = [(spec_x[i], average_spectrum[i]) for i in xrange(0, len_spectrum)]

max_spectra = sorted(average_spectrum_pairs, key=itemgetter(1), reverse=True)[:50]
local_max_spectra = argrelextrema(average_spectrum, np.greater)[0]
max_spectra_toplot = []
for i in xrange(0, len(max_spectra)):
	if (max_spectra[i][0]-1) in local_max_spectra:
		max_spectra_toplot.append(max_spectra[i])
max_spectra_toplot = sorted(max_spectra_toplot, key=itemgetter(0))

xtick_spread = len_spectrum / 100
max_spectra_ticks = [max_spectra_toplot[0]]
for i in xrange(1, len(max_spectra_toplot)):
	if max_spectra_toplot[i][0] - max_spectra_ticks[len(max_spectra_ticks)-1][0] < xtick_spread:
		if max_spectra_toplot[i][1] > max_spectra_ticks[len(max_spectra_ticks)-1][1]:
			max_spectra_ticks[len(max_spectra_ticks)-1] = max_spectra_toplot[i]
	else:
		max_spectra_ticks.append(max_spectra_toplot[i])

max_x = [max_spectra_toplot[i][0] for i in xrange(0, len(max_spectra_toplot))]
max_y = [max_spectra_toplot[i][1] for i in xrange(0, len(max_spectra_toplot))]
max_x_ticks = [max_spectra_ticks[i][0] for i in xrange(0, len(max_spectra_ticks))]
max_y_ticks = [max_spectra_ticks[i][1] for i in xrange(0, len(max_spectra_ticks))]
if max_x_ticks[0] > xtick_spread:
	max_x_ticks = [1] + max_x_ticks
	max_y_ticks = [average_spectrum[0]] + max_y_ticks
if max_x_ticks[len(max_x_ticks)-1] < len_spectrum - xtick_spread + 1:
	max_x_ticks = max_x_ticks + [len_spectrum+1]
	max_y_ticks = max_y_ticks + [average_spectrum[len_spectrum-1]]


fig = plt.figure(figsize=figsize_spectrum)
ax = fig.add_subplot(1,1,1)
ax.set_xlim(1, len_spectrum+1)
ax.plot(spec_x, average_spectrum, '-', color="blue")
ax.plot(max_x, max_y, 'bo', color="red")
ax.vlines(max_x, 0, max_y, color="black", linestyle="dotted")
ax.set_xticks(max_x_ticks)
ax.set_xticklabels(max_x_ticks, rotation='vertical')
fig.savefig('reports/latex/pics/mean_spectrum.pdf', format='pdf', bbox_inches='tight')
plt.close()

### read eigenvalues, eigenvectors, and coords
# eigvals_raw = genfromtxt('data/RB3.mat.val.csv', delimiter=';').astype(float)
# eigvecs_raw = genfromtxt('data/RB3.mat.vec.csv', delimiter=';').astype(float)
# coords_raw = genfromtxt('data/RB3.mat.coords.csv', delimiter=';').astype(float)
eigvals_raw = genfromtxt(args.e[0], delimiter=';').astype(float)
eigvecs_raw = genfromtxt(args.v[0], delimiter=';').astype(float)
coords_raw = genfromtxt(args.c[0], delimiter=';').astype(float)
eigvecs_raw[:,0] = eigvecs_raw[:,0].astype(int)
eigvecs_raw[:,1] = eigvecs_raw[:,1].astype(int)
coords_raw[:,0] = coords_raw[:,0] - min(coords_raw[:,0])
coords_raw[:,1] = coords_raw[:,1] - min(coords_raw[:,1])

rawvec_length = int(max(eigvecs_raw[:,1]))
vec_length = coords_raw.shape[0]
num_pixels = len(coords_raw)

# # sum of ions
sum_ions = np.sum(mat[spec_varname], axis=0)
fig = plt.figure()
image = np.zeros(shape = (int(max(coords_raw[:,0])+1), int(max(coords_raw[:,1])+1)))
colors = [('white')] + [(cm.jet(ind_color)) for ind_color in xrange(1,256)]
new_map = matplotlib.colors.LinearSegmentedColormap.from_list('new_map', colors, N=256)
for i in xrange(0, coords_raw.shape[0]):
	image[ int(coords_raw[i,0]), int(coords_raw[i,1]) ] = sum_ions[i]
plt.imshow(image, cmap=new_map, interpolation="None")
ax2 = fig.add_axes([0.92, 0.2, 0.03, 0.6])
cb = ColorbarBase(ax2, cmap=new_map, spacing='proportional', boundaries=linspace(0, max(sum_ions)), format='%1i')
fig.savefig('reports/latex/pics/sum_ions.pdf', format='pdf', bbox_inches='tight')
plt.close()

### select first num_eigens eigvecs
start = 0
if (abs(eigvals_raw[0,0]) < 1e-05):
	start = 1
eigvals = eigvals_raw[start:(start+num_eigens), 0]
eigvecs = np.zeros(shape = (vec_length, num_eigens-start))
for i in xrange(0, num_eigens-start):
	eigvecs[:, i] = eigvecs_raw[(i+start)*rawvec_length:(i+start)*rawvec_length + vec_length, 2]

white_eigvecs = vq.whiten(eigvecs)

for cur_num_centers in xrange(2, num_centers+1):
	### k-means clustering
	(centers, distortion) = vq.kmeans(white_eigvecs, cur_num_centers)
	labels = vq.vq(white_eigvecs, centers)[0]+1
	min_label = min(labels)
	cur_num_centers_res = max(labels)

	## color map for segmentation
	cmap = get_colormap(cur_num_centers)
	bounds = linspace(0,cur_num_centers_res+1,cur_num_centers_res+2)
	norm = BoundaryNorm(bounds, cmap.N)
	mappable = cm.ScalarMappable(norm=norm, cmap=cmap)
	## and a color map for heatmaps
	colors = [('white')] + [(cm.jet(ind_color)) for ind_color in xrange(1,256)]
	new_map = matplotlib.colors.LinearSegmentedColormap.from_list('new_map', colors, N=256)

	## total intensity by label
	total_intensities = np.zeros([len_spectrum, cur_num_centers])
	for l in xrange(0, num_pixels):
		total_intensities[:, labels[l]-1] = total_intensities[:, labels[l]-1] + mat[spec_varname][:, l]
	for i in xrange(0, cur_num_centers_res):
		total_intensities[:, i] = total_intensities[:, i] / sum(labels == i+1)
	spec_colors = np.argmax(total_intensities, axis=1)+1

	## total intensity in pixel by cluster
	total_bycluster = np.zeros([num_pixels, cur_num_centers])
	for i in xrange(0, cur_num_centers):
		total_bycluster[:, i] = np.sum( mat[spec_varname][spec_colors == i+1, :], axis=0 ) / sum(spec_colors == i+1)

	## intensities
	image_bycluster = np.zeros(shape = (cur_num_centers, cur_num_centers))
	for i in xrange(0, cur_num_centers):
		for j in xrange(0, cur_num_centers):
			image_bycluster[i, j] = np.mean( mat[spec_varname][spec_colors == i+1, :][:, labels == j+1] )

	fig = plt.figure(figsize=(10,10))
	# plt.tick_params(axis='both', which='major', labelsize="larger")
	aximage = plt.imshow(image_bycluster, cmap=cm.jet, interpolation="None", )
	plt.text(-0.75, -0.3, "pixels", horizontalalignment='center', verticalalignment='center', color="black", fontname="Sans", weight="bold")
	plt.text(cur_num_centers-.75, cur_num_centers-.35, "m/z values", horizontalalignment='center', verticalalignment='center', color="black", fontname="Sans", weight="bold")
	for i in xrange(0, cur_num_centers):
		# plt.text(i, -.75, ("%d\n%.2f%%" % ( sum(spec_colors == i+1), (100*sum(spec_colors == i+1) / float(len_spectrum)) )), horizontalalignment='center', verticalalignment='center', color="black", fontname="Sans", weight="bold")
		# plt.text(cur_num_centers-.2, i, ("%d\n%.2f%%" % ( sum(labels == i+1), (100*sum(labels == i+1) / float(num_pixels)) )), horizontalalignment='center', verticalalignment='center', color="black", fontname="Sans", weight="bold")
		for j in xrange(0, cur_num_centers):
			plt.text(j, i, ("%.02f" % image_bycluster[i,j]), horizontalalignment='center', verticalalignment='center', color="black", size="larger", fontname="Sans", weight="bold",
				path_effects=[pe.withStroke(linewidth=3, foreground="w")])
	ax2 = fig.add_axes([0.99, 0.2, 0.03, 0.6])
	Colorbar(ax2, aximage)
	fig.savefig('reports/latex/pics/' + str(cur_num_centers) + '_cluster_mean_intensities.pdf', format='pdf', bbox_inches='tight')
	plt.close()

	## numbers
	image_bycluster = np.zeros(shape = (cur_num_centers, cur_num_centers))
	for i in xrange(0, cur_num_centers):
		for j in xrange(0, cur_num_centers):
			shape = mat[spec_varname][spec_colors == i+1, :][:, labels == j+1].shape
			image_bycluster[i, j] = shape[0] * shape[1]

	fig = plt.figure(figsize=(12,10))
	aximage = plt.imshow(image_bycluster, cmap=cm.jet, interpolation="None", norm=matplotlib.colors.LogNorm())
	plt.text(-0.75, -0.3, "pixels", horizontalalignment='center', verticalalignment='center', color="black", fontname="Sans", weight="bold")
	plt.text(cur_num_centers-.75, cur_num_centers-.35, "m/z values", horizontalalignment='center', verticalalignment='center', color="black", fontname="Sans", weight="bold")
	for i in xrange(0, cur_num_centers):
		# plt.text(i, -.75, ("%d\n%.2f%%" % ( sum(spec_colors == i+1), (100*sum(spec_colors == i+1) / float(len_spectrum)) )), horizontalalignment='center', verticalalignment='center', color="black", fontname="Sans", weight="bold")
		# plt.text(cur_num_centers-.2, i, ("%d\n%.2f%%" % ( sum(labels == i+1), (100*sum(labels == i+1) / float(num_pixels)) )), horizontalalignment='center', verticalalignment='center', color="black", fontname="Sans", weight="bold")
		for j in xrange(0, cur_num_centers):
			plt.text(j, i, ("{0:,}\n{1:.2f}%".format(int(image_bycluster[i,j]), 100*image_bycluster[i,j] / float(num_pixels*len_spectrum))), horizontalalignment='center', verticalalignment='center', color="black", size="larger", fontname="Sans", weight="bold",
				path_effects=[pe.withStroke(linewidth=3, foreground="w")])
	ax2 = fig.add_axes([0.92, 0.2, 0.03, 0.6])
	Colorbar(ax2, aximage)
	fig.savefig('reports/latex/pics/' + str(cur_num_centers) + '_cluster_num_values.pdf', format='pdf', bbox_inches='tight')
	plt.close()

	## matrix view of the dataset
	mz_ordering = sorted(range(0, len_spectrum), key=lambda x: (spec_colors[x], -total_intensities[x, spec_colors[x]-1]) )
	pixel_ordering = sorted(range(0, num_pixels), key=lambda x: (labels[x], -total_bycluster[x, labels[x]-1] ) )
	image = np.zeros(shape = (num_pixels, len_spectrum))
	for i in xrange(0, len_spectrum):
		image[:, i] =  mat[spec_varname][mz_ordering[i], pixel_ordering]

	# fig = plt.figure(figsize=(10,10))
	# aximage = plt.imshow(image, cmap=cm.jet, interpolation="None")
	# # props = dict(boxstyle='round', facecolor=my_palette_array[i])
	# # plt.text(5, 5, str(i), horizontalalignment='left', verticalalignment='top', backgroundcolor=my_palette_array[i], color=my_palette_array[i], size="x-large", fontname="Sans", weight="bold",
	# # 	path_effects=[pe.withStroke(linewidth=3, foreground="w")], bbox=props)
	# ax2 = fig.add_axes([0.92, 0.2, 0.03, 0.6])
	# Colorbar(ax2, aximage)
	# fig.savefig('reports/latex/pics/' + str(cur_num_centers) + '_matrix_view.pdf', format='pdf', bbox_inches='tight')
	# plt.close()

	## average spectrum by cluster
	for i in xrange(1, cur_num_centers_res+1):
		fig = plt.figure(figsize=figsize_spectrum)
		ax = fig.add_subplot(1,1,1)
		ax.set_xlim(1, len_spectrum+1)
		ax.set_ylim(0, max(total_intensities[:, i-1])+1)
		# ax.plot(spec_x, average_spectrum, '-', color="gray")
		ax.plot(spec_x, total_intensities[:, i-1], linewidth=1.5, color=my_palette_array[i])
		props = dict(boxstyle='round', facecolor=my_palette_array[i])
		plt.text(.02*len_spectrum, .97*(max(total_intensities[:, i-1])+1), str(i),
			horizontalalignment='left', verticalalignment='top', backgroundcolor=my_palette_array[i], color=my_palette_array[i], size="x-large", fontname="Sans", weight="bold",
			path_effects=[pe.withStroke(linewidth=3, foreground="w")], bbox=props)
		fig.savefig('reports/latex/pics/' + str(cur_num_centers) + '_' + str(i) + '_mean_spectrum.pdf', format='pdf', bbox_inches='tight')
		plt.close()

	## colored spectrum picture
	fig = plt.figure(figsize=figsize_spectrum)
	ax = fig.add_subplot(1,1,1)
	ax.set_xlim(1, len_spectrum+1)
	ax.set_ylim(-1, max(average_spectrum))
	ax.plot(spec_x, average_spectrum, '-', color="gray")
	for i in xrange(1, cur_num_centers_res+1):
		ax.scatter([spec_x[j] for j in np.where(spec_colors == i)[0]], [average_spectrum[j] for j in np.where(spec_colors == i)[0]], color=my_palette_array[i])
		ax.scatter([spec_x[j] for j in np.where(spec_colors == i)[0]], [-0.5 for j in np.where(spec_colors == i)[0]], color=my_palette_array[i])
		# ax.plot([spec_x[j] for j in np.where(spec_colors == i)[0]], [average_spectrum[j] for j in np.where(spec_colors == i)[0]], 'o', color=my_palette_array[i], linewidth='0')
		# ax.plot([spec_x[j] for j in np.where(spec_colors == i)[0]], [0 for j in np.where(spec_colors == i)[0]], 's', color=my_palette_array[i], linewidth='0')
	ax2 = fig.add_axes([0.92, 0.1, 0.02, 0.8])
	cb = ColorbarBase(ax2, cmap=cmap, norm=norm, spacing='proportional', ticks=bounds, boundaries=bounds, format='%1i')
	fig.savefig('reports/latex/pics/' + str(cur_num_centers) + '_mean_spectrum.pdf', format='pdf', bbox_inches='tight')
	plt.close()

	for i in xrange(1, cur_num_centers+1):
		fig = plt.figure()
		image = np.zeros(shape = ( int(max(coords_raw[:,0])+1), int(max(coords_raw[:,1])+1) ) )
		for j in xrange(0, num_pixels):
			image[ int(coords_raw[j,0]), int(coords_raw[j,1]) ] = total_bycluster[j, i-1]
		aximage = plt.imshow(image, cmap=new_map, interpolation="None")
		props = dict(boxstyle='round', facecolor=my_palette_array[i])
		plt.text(5, 5, str(i), horizontalalignment='left', verticalalignment='top', backgroundcolor=my_palette_array[i], color=my_palette_array[i], size="x-large", fontname="Sans", weight="bold",
			path_effects=[pe.withStroke(linewidth=3, foreground="w")], bbox=props)
		ax2 = fig.add_axes([0.92, 0.2, 0.03, 0.6])
		Colorbar(ax2, aximage)
		fig.savefig('reports/latex/pics/' + str(cur_num_centers) + '_' + str(i) + '_cluster_heatmap.pdf', format='pdf', bbox_inches='tight')
		plt.close()

	## segmentation map picture
	fig = plt.figure()
	image = np.zeros(shape = ( int(max(coords_raw[:,0]))+1, int(max(coords_raw[:,1]))+1))
	for i in xrange(0, coords_raw.shape[0]):
		image[ int(coords_raw[i,0]), int(coords_raw[i,1]) ] = labels[i]
	plt.imshow(image, cmap=cmap, norm=norm, interpolation="None")
	ax2 = fig.add_axes([0.92, 0.2, 0.03, 0.6])
	cb = ColorbarBase(ax2, cmap=cmap, norm=norm, spacing='proportional', ticks=bounds, boundaries=bounds, format='%1i')
	fig.savefig('reports/latex/pics/' + str(cur_num_centers) + '_segmentation.pdf', format='pdf', bbox_inches='tight')
	plt.close()


## LaTeX template
template = file('reports/latex/template.tex', 'r').read()
res = template % {'report_name' : args.t[0].replace('_', "\_"), 'num_pixels' : num_pixels, 'xdim' : image.shape[0], 'ydim' : image.shape[1], 'len_spectrum' : len_spectrum, 'num_clusters' : num_centers }
file('reports/latex/tmp.tex', 'w').write(res)
