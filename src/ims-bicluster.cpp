#include <math.h>
#include <omp.h>
#include <getopt.h>
#include <math.h>
#include <time.h>

#include <cstdlib>
#include <functional>
#include <iostream>
#include <fstream>
#include <sstream>
#include <algorithm>
#include <string>
#include <iomanip>
#include <sstream>
#include <vector>
#include <queue>
#include <list>
#include <ctime>

#include <boost/bind.hpp>
#include <boost/filesystem.hpp>
#include <boost/filesystem/fstream.hpp>
#include <boost/unordered_map.hpp>
#include <boost/algorithm/string.hpp>
#include <boost/lexical_cast.hpp>

#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_eigen.h>
#include <omp.h>

#include <matio.h>

#define no_argument 0
#define required_argument 1 
#define optional_argument 2

#include "logging.hpp"

using namespace std;

string input_filename = "";
bool matlab_input = false;

// weight of spatial edges in the incidence matrix
double alpha = 0.5;
// edge weight between neighboring points
double max_force = 120;
// distance inside which we introduce spatial edges
double threshold = 1.5;

int help(char *argv[]) {
    cout << "Usage: " << argv[0] << " [--mat] --input=input_file" << endl;
    return 0;
}

int main(int argc, char *argv[]) {
    const struct option longopts[] = {
        {"help",            no_argument,        0, 'h'},
        {"input",           required_argument,  0, 'i'},
        {"mat",             no_argument,        0, 'm'},
        {"alpha",           required_argument,  0, 'a'},
        {0,0,0,0},
    };

    if (argc < 1) return help(argv);
    timer.reset();
    int index, iarg=0;
    opterr=1;    
    while (iarg != -1) {
        iarg = getopt_long(argc, argv, "fvh", longopts, &index);
        switch (iarg) {
            case 'h':   return help(argv);              break;
            case 'i':   input_filename = optarg;        break;
            case 'a':   alpha = atof(optarg);           break;
            case 'm':   matlab_input = true;            break;
        }
    }

    // weight of bipartite edges in the incidence matrix
    double beta = 1 - alpha;

    mat_t *matfp;
    matvar_t *matvar;

    uint num_pixels, len_spectrum;
    double **spectra;
    double *maxima, *specsdiag, *pixdiag;
    
    if (matlab_input) {
        LOG("Reading input from " << input_filename << " as a Matlab file...");
        matfp = Mat_Open(input_filename.c_str(), MAT_ACC_RDONLY);

        if ( NULL == matfp ) {
            fprintf(stderr,"Error creating MAT file \"matfile5.mat\"!\n");
            return EXIT_FAILURE;
        }

        matvar = Mat_VarRead(matfp, "spectra");
        num_pixels = matvar->dims[0];
        len_spectrum = matvar->dims[1];

        spectra = allocate_2d_with_default<double>(num_pixels, len_spectrum, 0);
        maxima = allocate_1d_with_default<double>(num_pixels, 0);
        specsdiag = allocate_1d_with_default<double>(len_spectrum, 0);
        pixdiag = allocate_1d_with_default<double>(num_pixels, 0);

        double val;
        LOG("\t\t" << ((double*)(matvar->data))[0] << " " << ((double*)(matvar->data))[1]);
        for (uint i=0; i<num_pixels; ++i) {
            for (uint j=0; j<len_spectrum; ++j) {
                val = ((double*)(matvar->data))[num_pixels*j+i];
                spectra[i][j] = val;
                specsdiag[j] += val;
                pixdiag[i] += val;
                if (val > maxima[i]) maxima[i] = val;
            }
        }
        Mat_VarFree(matvar);
        Mat_Close(matfp);
    }

    uint k = (uint)(8*num_pixels/(double)9);
    double r = 0.20;

    // find maximal elements
    LOG("Finding maximal elements...");
    uint **U = allocate_2d_with_default<uint>(num_pixels, len_spectrum, 0);
    vector<uint> indices(len_spectrum);
    for (uint j=0; j<len_spectrum; ++j) { indices.push_back(j); }
    for (uint i=0; i<num_pixels; ++i) {
        sort(begin(indices), end(indices), [&](const uint & j1, const uint & j2) {
            return spectra[i][j1] > spectra[i][j2];
        });
        for (uint j=0; j<k; ++j) {
            U[i][indices[j]] = 1;
        }
        // LOG("\t" << spectra[i][indices[0]] << "\t" << spectra[i][indices[1]] << "\t" << spectra[i][indices[2]]);
    }

    LOG("Filling incidence matrix...");
    double **W = allocate_2d_with_default<double>(num_pixels, num_pixels, 0);
    double *on_diag = allocate_1d_with_default<double>(num_pixels, 0);
    double *little = allocate_1d_with_default<double>(num_pixels, 0);
    

    delete maxima;
    delete specsdiag;
    delete pixdiag;
    delete_2d<double>(num_pixels, spectra);

    LOG("All done.");
    return 0;
}


    // double data[] = { 1.0  , 1/2.0, 1/3.0, 1/4.0,
    //                 1/2.0, 1/3.0, 1/4.0, 1/5.0,
    //                 1/3.0, 1/4.0, 1/5.0, 1/6.0,
    //                 1/4.0, 1/5.0, 1/6.0, 1/7.0 };

    // gsl_matrix_view m = gsl_matrix_view_array (data, 4, 4);

    // gsl_vector *eval = gsl_vector_alloc (4);
    // gsl_matrix *evec = gsl_matrix_alloc (4, 4);

    // gsl_eigen_symmv_workspace * w =  gsl_eigen_symmv_alloc (4);
    // gsl_eigen_symmv (&m.matrix, eval, evec, w);
    // gsl_eigen_symmv_free (w);
    // gsl_eigen_symmv_sort (eval, evec, GSL_EIGEN_SORT_ABS_ASC);

    // LOG("Sample output from gsl_eigen:")  
    // for (uint i = 0; i < 4; ++i) {
    //     double eval_i = gsl_vector_get (eval, i);
    //     gsl_vector_view evec_i = gsl_matrix_column (evec, i);

    //     LOG ("\teigenvalue = " << eval_i);
    //     LOG ("\teigenvector = ");
    //     gsl_vector_fprintf (stdout, &evec_i.vector, "%g");
    // }

