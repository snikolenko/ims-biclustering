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

#include <matio.h>

#define no_argument 0
#define required_argument 1 
#define optional_argument 2

#include "logging.hpp"

#include "arpack++/areig.h"
// #include "lnmatrxc.h"

using namespace std;

typedef double t_ims_real;

string input_filename = "";
bool matlab_input = false;
bool matlab_input_2 = false;

// weight of spatial edges in the incidence matrix
t_ims_real alpha = 0.5;
// edge weight between neighboring points
t_ims_real max_force = 120;
// distance inside which we introduce spatial edges
t_ims_real threshold = 1.5 * 1.5;
// denominator threshold
t_ims_real r = 0.20;
// how many eigenvalues to find
int num_eigens = 10;

int help(char *argv[]) {
    cout << "Usage: " << argv[0] << " [--mat] --input=input_file" << endl;
    return 0;
}

int main(int argc, char *argv[]) {
    const struct option longopts[] = {
        {"help",            no_argument,        0, 'h'},
        {"input",           required_argument,  0, 'i'},
        {"mat",             no_argument,        0, 'm'},
        {"mat2",            no_argument,        0, 'n'},
        {"alpha",           required_argument,  0, 'a'},
        {"eigens",          required_argument,  0, 'e'},
        {"r",               required_argument,  0, 'r'},
        {"maxforce",        required_argument,  0, 'f'},
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
            case 'r':   r = atof(optarg);               break;
            case 'f':   max_force = atof(optarg);       break;
            case 'e':   num_eigens = atoi(optarg);      break;
            case 'm':   matlab_input = true;            break;
            case 'n':   matlab_input_2 = true;          break;
        }
    }


    // weight of bipartite edges in the incidence matrix
    t_ims_real beta = 1 - alpha;

    mat_t *matfp;
    matvar_t *matvar;

    uint len_spectrum, num_pixels;
    t_ims_real **spectra;
    t_ims_real *maxima, *specsdiag, *pixdiag;
    t_ims_real *xcoord, *ycoord;
    
    if (matlab_input || matlab_input_2) {
        LOG("Reading input from " << input_filename << " as a Matlab file...");
        matfp = Mat_Open(input_filename.c_str(), MAT_ACC_RDONLY);

        if ( NULL == matfp ) {
            fprintf(stderr,"Error creating MAT file \"matfile5.mat\"!\n");
            return EXIT_FAILURE;
        }

        if (matlab_input) {
            matvar = Mat_VarRead(matfp, "x_printed");
            num_pixels = matvar->dims[0];
            xcoord = allocate_1d_with_default<t_ims_real>(num_pixels, 0);
            for (uint j=0; j<num_pixels; ++j) {
                xcoord[j] = ((t_ims_real*)(matvar->data))[j];
            }

            matvar = Mat_VarRead(matfp, "y_printed");
            ycoord = allocate_1d_with_default<t_ims_real>(num_pixels, 0);
            for (uint j=0; j<num_pixels; ++j) {
                ycoord[j] = ((t_ims_real*)(matvar->data))[j];
            }
            
            matvar = Mat_VarRead(matfp, "spectra");
            len_spectrum = matvar->dims[0];
            LOG("Read " << matvar->dims[0] << " x " << matvar->dims[1] << " matrix.");
        }

        if (matlab_input_2) {
            matvar = Mat_VarRead(matfp, "coords");
            char * const * c = Mat_VarGetStructFieldnames(matvar);
            matvar_t * matvar_x = Mat_VarGetStructFieldByName(matvar, "x", 0);
            num_pixels = matvar_x->dims[1];
            xcoord = allocate_1d_with_default<t_ims_real>(num_pixels, 0);
            for (uint j=0; j<num_pixels; ++j) {
                xcoord[j] = ((t_ims_real*)(matvar_x->data))[j];
            }
            matvar_t * matvar_y = Mat_VarGetStructFieldByName(matvar, "y", 0);
            ycoord = allocate_1d_with_default<t_ims_real>(num_pixels, 0);
            for (uint j=0; j<num_pixels; ++j) {
                ycoord[j] = ((t_ims_real*)(matvar_y->data))[j];
            }
            
            matvar = Mat_VarRead(matfp, "SP");
            len_spectrum = matvar->dims[0];
            LOG("Read " << matvar->dims[0] << " x " << matvar->dims[1] << " matrix.");
        }
    LOG("Printing pixel coords to " << input_filename << ".coords.csv");
    ofstream ofs_coord(input_filename + ".coords.csv");
    for (uint i=0; i<num_pixels; ++i) {
        ofs_coord << xcoord[i] << ";" << ycoord[i] << "\n";
    }
    ofs_coord.close();


        spectra = allocate_2d_with_default<t_ims_real>(num_pixels, len_spectrum, 0);
        maxima = allocate_1d_with_default<t_ims_real>(num_pixels, 0);
        specsdiag = allocate_1d_with_default<t_ims_real>(len_spectrum, 0);
        pixdiag = allocate_1d_with_default<t_ims_real>(num_pixels, 0);

        t_ims_real val;
        for (uint i=0; i<len_spectrum; ++i) {
            for (uint j=0; j<num_pixels; ++j) {
                val = ((t_ims_real*)(matvar->data))[len_spectrum*j+i];
                spectra[j][i] = val;
                specsdiag[i] += val;
                pixdiag[j] += val;
                if (val > maxima[j]) maxima[j] = val;
            }
        }
        Mat_VarFree(matvar);
        Mat_Close(matfp);
    }

    uint k = (uint)(8*len_spectrum/(t_ims_real)9);

    // find maximal elements
    LOG("Finding maximal elements...");
    uint **U = allocate_2d_with_default<uint>(len_spectrum, num_pixels, 0);
    vector<uint> indices(num_pixels);
    for (uint j=0; j<num_pixels; ++j) { indices.push_back(j); }
    for (uint i=0; i<len_spectrum; ++i) {
        sort(begin(indices), end(indices), [&](const uint & j1, const uint & j2) {
            return spectra[i][j1] > spectra[i][j2];
        });
        for (uint j=0; j<k; ++j) {
            U[i][indices[j]] = 1;
        }
        // LOG("\t" << spectra[i][indices[0]] << "\t" << spectra[i][indices[1]] << "\t" << spectra[i][indices[2]]);
    }

    LOG("Filling incidence matrix...");
    t_ims_real **W = allocate_2d_with_default<t_ims_real>(num_pixels, num_pixels, 0);
    t_ims_real *on_diag = allocate_1d_with_default<t_ims_real>(num_pixels, 0);
    t_ims_real *little = allocate_1d_with_default<t_ims_real>(num_pixels, 0);

    for (uint j1=0; j1<num_pixels; ++j1) {
        for (uint j2=0; j2<j1; ++j2) {
            t_ims_real dist = (xcoord[j1]-xcoord[j2])*(xcoord[j1]-xcoord[j2]) + (ycoord[j1]-ycoord[j2])*(ycoord[j1]-ycoord[j2]);
            if (dist <= threshold) {
                uint common_pixels = 0;
                for (uint i=0; i<len_spectrum; ++i) {
                    if (U[i][j1] == 1 && U[i][j2] == 1) {
                        ++common_pixels;
                    }
                }
                W[j1][j2] = alpha * max_force * common_pixels / ( (len_spectrum - k) * dist );
                W[j2][j1] = W[j1][j2];
                on_diag[j1] += W[j1][j2];
                on_diag[j2] += W[j1][j2];
                if (common_pixels <= r * (len_spectrum - k)) {
                    little[j1] += W[j1][j2];
                    little[j2] += W[j1][j2];
                }
            }
        }
    }

    LOG("Allocating big matrix...");
    // t_ims_real **bigL = allocate_2d_with_default<t_ims_real>(num_pixels + len_spectrum, num_pixels + len_spectrum, 0);
    // t_ims_real **bigW = allocate_2d_with_default<t_ims_real>(num_pixels + len_spectrum, num_pixels + len_spectrum, 0);
    // arma::mat bigL = arma::zeros<arma::mat>(num_pixels + len_spectrum, num_pixels + len_spectrum);
    // Eigen::MatrixXf bigL(num_pixels + len_spectrum, num_pixels + len_spectrum);
    // gsl_matrix_set_zero(bigL);
    // gsl_matrix *bigW = gsl_matrix_alloc (num_pixels + len_spectrum, num_pixels + len_spectrum);
    // gsl_matrix_set_zero(bigW);

    int n = num_pixels + len_spectrum;
    int nnz = num_pixels * (num_pixels + len_spectrum) + num_pixels * len_spectrum + len_spectrum;
    int *irow = new int[nnz];
    int *pcol = new int[n+1];
    t_ims_real *A = new t_ims_real[nnz];

    LOG("Filling big matrix in ARPACK's column major format...");

    uint m=0;
    pcol[0] = 0;
    t_ims_real w_inv;
    for (uint j=0; j<n; ++j) { // column by column
        w_inv = 1.0 / ( (j < num_pixels) ? (little[j] + beta * pixdiag[j]) : (beta * specsdiag[j-num_pixels]) ); // will divide by W
        if (j < num_pixels) {
            for (uint i=0; i<num_pixels; ++i) { // top left block W
                irow[m] = i;
                A[m++] = ( (i == j ? (on_diag[j] + beta * pixdiag[j]) : 0) - W[i][j] ) * w_inv;
            }
            for (uint i=0; i<len_spectrum; ++i) { // bottom left block beta*spectra
                irow[m] = num_pixels + i;
                A[m++] = -beta * spectra[j][i] * w_inv;
            }
        } else {
            for (uint i=0; i<num_pixels; ++i) { // top right block beta*transpose(spectra)
                irow[m] = i;
                A[m++] = -beta * spectra[i][j-num_pixels] * w_inv;
            }
            // bottom right block -- diagonal
            irow[m] = j;
            A[m++] = beta * specsdiag[j-num_pixels] * w_inv;
        }
        pcol[j+1] = m;
    }

    LOG("Solving for " << num_eigens << " eigenvalues with ARPACK...");

    int nconv;
    t_ims_real *EigValR = new t_ims_real[num_eigens];
    t_ims_real *EigValI = new t_ims_real[num_eigens];
    t_ims_real *EigVec = new t_ims_real[n * num_eigens];

    nconv = AREig(EigValR, EigValI, EigVec, n, nnz, A, irow, pcol, num_eigens, "SR");
    LOG("Eigenvalues:");
    for (uint i=0; i<nconv; ++i) {
        LOG("\tlambda[" << (i+1) << "] = " << EigValR[i] << (EigValI[i] >= 0.0 ? '+' : '-') << fabs(EigValI[i]));
    }

    LOG("Printing eigenvalues to " << input_filename << ".val.csv");
    ofstream ofs_val(input_filename + ".val.csv");
    for (uint i=0; i<nconv; ++i) {
        ofs_val << EigValR[i] << ";" << EigValI[i] << "\n";
    }
    ofs_val.close();

    LOG("Printing eigenvectors to " << input_filename << ".vec.csv");
    ofstream ofs_vec(input_filename + ".vec.csv");
    for (uint i=0; i<nconv; ++i) {
        for (uint j=0; j<n; ++j) {
            ofs_vec << (i+1) << ";" << (j+1) << ";" << EigVec[i*n+j] << "\n";
        }
    }
    ofs_vec.close();

    LOG("Printing pixel coords to " << input_filename << ".coords.csv");
    ofstream ofs_coord(input_filename + ".coords.csv");
    for (uint i=0; i<num_pixels; ++i) {
        ofs_coord << xcoord[i] << ";" << ycoord[i] << "\n";
    }
    ofs_coord.close();

    delete maxima;
    delete specsdiag;
    delete pixdiag;
    delete on_diag;
    delete little;
    delete [] EigValR;
    delete [] EigValI;
    delete [] EigVec;
    delete [] irow;
    delete [] pcol;
    delete [] A;
    delete_2d<t_ims_real>(len_spectrum, spectra);
    delete_2d<t_ims_real>(num_pixels, W);

    LOG("All done.");
    return 0;
}
