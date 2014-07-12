import matplotlib.pyplot as plt
from scipy.cluster.vq import kmeans2
from scipy.cluster.vq import kmeans
import sklearn.cluster as cluster

num_centers = 10

### read eigenvalues and eigenvectors
eigvals_raw = genfromtxt('data/RB3.mat.val.csv', delimiter=';').astype(float)
eigvecs_raw = genfromtxt('data/RB3.mat.vec.csv', delimiter=';').astype(float)
eigvecs_raw[:,0] = eigvecs_raw[:,0].astype(int)
eigvecs_raw[:,1] = eigvecs_raw[:,1].astype(int)

### read x and y coords
coords_raw = genfromtxt('data/RB3.mat.coords.csv', delimiter=';').astype(float)
coords_raw[:,0] = coords_raw[:,0] - min(coords_raw[:,0])
coords_raw[:,1] = coords_raw[:,1] - min(coords_raw[:,1])

vec_length = max(eigvecs_raw[:,1])
num_pixels = len(coords_raw)

### select first num_centers eigvecs
start = 0
if (abs(eigvals_raw[0,0]) < 1e-05):
	start = 1
eigvals = eigvals_raw[start:(start+num_centers), 0]
eigvecs = numpy.zeros(shape = (vec_length, num_centers))
for i in xrange(0, num_centers):
	eigvecs[:, i] = eigvecs_raw[(i+start)*vec_length:(i+start+1)*vec_length, 2]

### k-means clustering
(centers, labels) = kmeans2(eigvecs, num_centers, minit = 'random')

### spectral
clust = cluster.SpectralClustering(num_centers)
labels_spec = clust.fit_predict(eigvecs)
min_label = min(labels)

image = numpy.zeros(shape = (max(coords_raw[:,0])+1, max(coords_raw[:,1])+1))
for i in xrange(0, coords_raw.shape[0]):
	image[ coords_raw[i,0], coords_raw[i,1] ] = labels[i] - min_label + 1



## color map magic
cmap = plt.cm.Accent
cmaplist = [cmap(i) for i in range(cmap.N)]
cmaplist[0] = (1,1,1,1.0)
cmap = cmap.from_list('Custom cmap', cmaplist, cmap.N)
bounds = np.linspace(0,num_centers,num_centers+1)
norm = mpl.colors.BoundaryNorm(bounds, cmap.N)

fig = plt.figure()
plt.imshow(image, cmap=cmap, norm=norm, interpolation="None")
ax2 = fig.add_axes([0.92, 0, 0.03, 1])
cb = mpl.colorbar.ColorbarBase(ax2, cmap=cmap, norm=norm, spacing='proportional', ticks=bounds, boundaries=bounds, format='%1i')
# ax2.set_ylabel("Colorbar for the clusters")
fig.savefig('fig.pdf', format='pdf')
#plt.show()


