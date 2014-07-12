import matplotlib.pyplot as plt
import scipy.cluster.vq as vq
import sklearn.cluster as cluster
import argparse
from numpy import genfromtxt,zeros,linspace
from matplotlib.colors import BoundaryNorm
from matplotlib.colorbar import ColorbarBase

parser = argparse.ArgumentParser(description='Cluster eigenvectors and generate a pretty report.')
parser.add_argument('-t', metavar='TITLE', type=str, nargs='+', help='project title')
parser.add_argument('-e', metavar='VALUES', type=str, nargs='+', help='file of eigenvalues')
parser.add_argument('-v', metavar='VECTORS', type=str, nargs='+', help='file of eigenvectors')
parser.add_argument('-c', metavar='COORDS', type=str, nargs='+', help='file of coordinates')
parser.add_argument('-N', metavar='NUM_CLUSTERS', type=int, nargs='+', help='number of clusters')
parser.add_argument('-M', metavar='NUM_EIGENS', type=int, nargs='+', help='number of eigenvalues to consider')

args = parser.parse_args()

num_eigens = args.M[0]
num_centers = args.N[0]

# num_eigens = 10
# num_centers = 15
# eigvals_raw = genfromtxt('data/RB3.mat.val.csv', delimiter=';').astype(float)
# eigvecs_raw = genfromtxt('data/RB3.mat.vec.csv', delimiter=';').astype(float)
# coords_raw = genfromtxt('data/RB3.mat.coords.csv', delimiter=';').astype(float)

### read eigenvalues, eigenvectors, and coords
eigvals_raw = genfromtxt(args.e[0], delimiter=';').astype(float)
eigvecs_raw = genfromtxt(args.v[0], delimiter=';').astype(float)
coords_raw = genfromtxt(args.c[0], delimiter=';').astype(float)
eigvecs_raw[:,0] = eigvecs_raw[:,0].astype(int)
eigvecs_raw[:,1] = eigvecs_raw[:,1].astype(int)
coords_raw[:,0] = coords_raw[:,0] - min(coords_raw[:,0])
coords_raw[:,1] = coords_raw[:,1] - min(coords_raw[:,1])

rawvec_length = max(eigvecs_raw[:,1])
vec_length = coords_raw.shape[0]
num_pixels = len(coords_raw)

### select first num_centers eigvecs
start = 0
if (abs(eigvals_raw[0,0]) < 1e-05):
	start = 1
eigvals = eigvals_raw[start:(start+num_eigens), 0]
eigvecs = zeros(shape = (vec_length, num_eigens))
for i in xrange(0, num_eigens):
	eigvecs[:, i] = eigvecs_raw[(i+start)*rawvec_length:(i+start)*rawvec_length + vec_length, 2]

white_eigvecs = vq.whiten(eigvecs)

### k-means clustering
(centers, distortion) = vq.kmeans(white_eigvecs, num_centers)
labels = vq.vq(white_eigvecs, centers)[0]

min_label = min(labels)
num_centers_res = max(labels)
image = zeros(shape = (max(coords_raw[:,0])+1, max(coords_raw[:,1])+1))
for i in xrange(0, coords_raw.shape[0]):
	image[ coords_raw[i,0], coords_raw[i,1] ] = labels[i] - min_label + 1

## color map magic
cmap = plt.cm.Accent
cmaplist = [cmap(i) for i in range(cmap.N)]
cmaplist[0] = (1,1,1,1.0)
cmap = cmap.from_list('Custom cmap', cmaplist, cmap.N)
bounds = linspace(0,num_centers_res,num_centers_res+1)
norm = BoundaryNorm(bounds, cmap.N)


## segmentation map picture
fig = plt.figure()
plt.imshow(image, cmap=cmap, norm=norm, interpolation="None")
ax2 = fig.add_axes([0.92, 0.2, 0.03, 0.6])
cb = ColorbarBase(ax2, cmap=cmap, norm=norm, spacing='proportional', ticks=bounds, boundaries=bounds, format='%1i')
fig.savefig('reports/latex/segmentation.pdf', format='pdf')


## LaTeX template
template = file('reports/latex/template.tex', 'r').read()
res = template % {'report_name' : args.t[0].replace('_',"\_")}
file('reports/latex/tmp.tex', 'w').write(res)
