ims-biclustering
================

Biclustering algorithms for processing imaging mass-spectrometry data.


=== Installation ===

To compute eigenvalues, IMS-Biclustering uses ARPACK; the ARPACK++ headers are included. The following libraries have to be installed:
	libblas libsuperlu libarpack libgfortran libgfortranbegin libnsl
To install them on an Ubuntu system, run
	sudo apt-get install libarpack2-dev libsuperlu3-dev gfortran libatlas libatlas-dev libblas-dev 

Reports generation uses LaTeX and python with numpy and matplotlib installed:
	sudo apt-get install latex texlive dvipng libfreetype6-dev libpng-dev
	sudo pip install matplotlib numpy scikit-learn statsmodels


=== Usage: bash scripts ===

There are two bash scripts for the IMS dataset processing cycle:

	./gen-all.sh fname
Expects to find input file in data/fname as a Matlab .mat file. Runs ims-bicluster, then runs gen-report to generate the report pdf.

	./gen-report.sh fname
Expects to find input file in data/fname as a Matlab .mat file together with the following files produced by ims-bicluster:
	fname.mat.coords.csv	-- coordinates of spectra on the image
	fname.mat.val.csv		-- eigenvalues
	fname.mat.vec.csv		-- eigenvectors
 Runs ims-bicluster, then runs gen-report to generate the report pdf.


=== Usage: ims-bicluster ===

	 ./bin/ims-bicluster --mat --eigens=n --input=ims-dataset.mat
Computes n largest eigenvalues for the matrix encoded as a Matlab .mat file with the following structure:
 * x_printed	-- vector of X coordinates of all pixels in the dataset (len=num_pixels)
 * y_printed	-- vector of Y coordinates of all pixels in the dataset (len=num_pixels)
 * spectra		-- matrix of intensities for each pixel (shape = num_pixels x len_spectrum)

	 ./bin/ims-bicluster --mat2 --eigens=n --input=ims-dataset.mat
Computes n largest eigenvalues for the matrix encoded as a Matlab .mat file with the following structure:
 * x	-- vector of X coordinates of all pixels in the dataset (len=num_pixels)
 * y	-- vector of Y coordinates of all pixels in the dataset (len=num_pixels)
 * SP	-- matrix of intensities for each pixel (shape = num_pixels x len_spectrum)


=== Usage: visualize.py ===

	python python/visualize.py -e data/fname.val.csv -v data/fname.vec.csv -c data/fname.coords.csv -N num_clusters -M num_eigenvalues
Performs clustering for all numbers of clusters from 2 to num_clusters and num_eigenvalues largest eigenvalues based on ims-bicluster processing results; generates a pdf report. The LaTex code of the last processed report is stored at reports/latex/tmp.tex.


