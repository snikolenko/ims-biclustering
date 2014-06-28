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

 #define ARMA_DONT_USE_CXX11
#include <armadillo>

#define no_argument 0
#define required_argument 1 
#define optional_argument 2

#include "logging.hpp"

#include <redsvd/redsvd.hpp>
#include <eigen3/Eigen/Dense>

using namespace std;

string input_filename = "";
bool matlab_input = false;

// weight of spatial edges in the incidence matrix
double alpha = 0.5;
// edge weight between neighboring points
double max_force = 120;
// distance inside which we introduce spatial edges
double threshold = 1.5 * 1.5;

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

    uint len_spectrum, num_pixels;
    double **spectra;
    double *maxima, *specsdiag, *pixdiag;
    double *xcoord, *ycoord;
    
    if (matlab_input) {
        LOG("Reading input from " << input_filename << " as a Matlab file...");
        matfp = Mat_Open(input_filename.c_str(), MAT_ACC_RDONLY);

        if ( NULL == matfp ) {
            fprintf(stderr,"Error creating MAT file \"matfile5.mat\"!\n");
            return EXIT_FAILURE;
        }

        matvar = Mat_VarRead(matfp, "x_printed");
        num_pixels = matvar->dims[0];
        xcoord = allocate_1d_with_default<double>(num_pixels, 0);
        for (uint j=0; j<num_pixels; ++j) {
            xcoord[j] = ((double*)(matvar->data))[j];
        }

        matvar = Mat_VarRead(matfp, "y_printed");
        ycoord = allocate_1d_with_default<double>(num_pixels, 0);
        for (uint j=0; j<num_pixels; ++j) {
            ycoord[j] = ((double*)(matvar->data))[j];
        }
        
        matvar = Mat_VarRead(matfp, "spectra");
        len_spectrum = matvar->dims[0];
        LOG(matvar->dims[0] << " " << matvar->dims[1]);

        spectra = allocate_2d_with_default<double>(num_pixels, len_spectrum, 0);
        maxima = allocate_1d_with_default<double>(num_pixels, 0);
        specsdiag = allocate_1d_with_default<double>(len_spectrum, 0);
        pixdiag = allocate_1d_with_default<double>(num_pixels, 0);

        double val;
        LOG("\t\t" << ((double*)(matvar->data))[0] << " " << ((double*)(matvar->data))[1]);
        for (uint i=0; i<len_spectrum; ++i) {
            for (uint j=0; j<num_pixels; ++j) {
                val = ((double*)(matvar->data))[len_spectrum*j+i];
                spectra[j][i] = val;
                specsdiag[i] += val;
                pixdiag[j] += val;
                if (val > maxima[j]) maxima[j] = val;
            }
        }
        Mat_VarFree(matvar);
        Mat_Close(matfp);
    }

    uint k = (uint)(8*len_spectrum/(double)9);
    double r = 0.20;

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
    double **W = allocate_2d_with_default<double>(num_pixels, num_pixels, 0);
    double *on_diag = allocate_1d_with_default<double>(num_pixels, 0);
    double *little = allocate_1d_with_default<double>(num_pixels, 0);

    for (uint j1=0; j1<num_pixels; ++j1) {
        for (uint j2=0; j2<j1; ++j2) {
            double dist = (xcoord[j1]-xcoord[j2])*(xcoord[j1]-xcoord[j2]) + (ycoord[j1]-ycoord[j2])*(ycoord[j1]-ycoord[j2]);
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

    LOG("Allocating big matrices...");
    // double **bigL = allocate_2d_with_default<double>(num_pixels + len_spectrum, num_pixels + len_spectrum, 0);
    // double **bigW = allocate_2d_with_default<double>(num_pixels + len_spectrum, num_pixels + len_spectrum, 0);
    arma::mat bigL = arma::zeros<arma::mat>(num_pixels + len_spectrum, num_pixels + len_spectrum);
    // Eigen::MatrixXf bigL(num_pixels + len_spectrum, num_pixels + len_spectrum);
    // gsl_matrix_set_zero(bigL);
    // gsl_matrix *bigW = gsl_matrix_alloc (num_pixels + len_spectrum, num_pixels + len_spectrum);
    // gsl_matrix_set_zero(bigW);

    LOG("Filling big matrices...");

    // top left block W
    for (uint j1=0; j1<num_pixels; ++j1) {
        for (uint j2=0; j2<num_pixels; ++j2) {
            // gsl_matrix_set(bigL, j1, j2, -W[j1][j2]);
            // bigL[j1][j2] = -W[j1][j2];
            // bigL.at(j1, j2) = -W[j1][j2];
            bigL(j1, j2) = -W[j1][j2];
        }
    }

    LOG("\ttop left done...");

    // top right block beta*transpose(spectra)
    for (uint j=0; j<num_pixels; ++j) {
        // LOG("\t\t" << j);
        for (uint i=0; i<len_spectrum; ++i) {
            // gsl_matrix_set(bigL, j, num_pixels + i, -beta * spectra[j][i]);
            // bigL[j][num_pixels + i] = -beta * spectra[j][i];
            // bigL.at(j, num_pixels + i) = -beta * spectra[j][i];
            bigL(j, num_pixels + i) = -beta * spectra[j][i];
        }
    }

    LOG("\ttop right done...");

    // bottom left block beta*spectra
    for (uint j=0; j<num_pixels; ++j) {
        for (uint i=0; i<len_spectrum; ++i) {
            // gsl_matrix_set(bigL, num_pixels + i, j, -beta * spectra[j][i]);
            // bigL[num_pixels + i][j] = -beta * spectra[j][i];
            bigL(num_pixels + i, j) = -beta * spectra[j][i];
        }
    }

    LOG("\tbottom left done...");

    // bottom right block remains zero?

    // now add bigD
    for (uint j=0; j<num_pixels; ++j) {
        // bigL[j][j] += on_diag[j] + beta * pixdiag[j];
        // bigW[j][j] += little[j] + beta * pixdiag[j];
        // gsl_matrix_set(bigL, j, j, gsl_matrix_get(bigL, j, j) + on_diag[j] + beta * pixdiag[j]);
        bigL(j, j) = bigL(j, j) + on_diag[j] + beta * pixdiag[j];
        // bigL.at(j, j) = bigL.at(j, j) + on_diag[j] + beta * pixdiag[j];
        // gsl_matrix_set(bigW, j, j, gsl_matrix_get(bigW, j, j) + little[j] + beta * pixdiag[j]);
    }

    for (uint i=0; i<len_spectrum; ++i) {
        // bigL[num_pixels + i][num_pixels + i] += beta * specsdiag[i];
        // bigW[num_pixels + i][num_pixels + i] += beta * specsdiag[i];
        // gsl_matrix_set(bigL, num_pixels + i, num_pixels + i, gsl_matrix_get(bigL, num_pixels + i, num_pixels + i) + beta * specsdiag[i]);
        // bigL.at(num_pixels + i, num_pixels + i) = bigL.at(num_pixels + i, num_pixels + i) + beta * specsdiag[i];
        bigL(num_pixels + i, num_pixels + i) = bigL(num_pixels + i, num_pixels + i) + beta * specsdiag[i];
        // gsl_matrix_set(bigW, num_pixels + i, num_pixels + i, gsl_matrix_get(bigW, num_pixels + i, num_pixels + i) + beta * specsdiag[i]);
    }

    LOG("\tdiagonal done...");

    // divide by bigW
    double w_inv;
    for (uint j=0; j<num_pixels+len_spectrum; ++j) {
        w_inv = 1.0 / ( (j < num_pixels) ? (little[j] + beta * pixdiag[j]) : (beta * specsdiag[j-num_pixels]) );
        for (uint i=0; i<num_pixels+len_spectrum; ++i) {
            bigL(i,j) = bigL(i,j) * w_inv;
        }
    }

    LOG("\tdivided by bigW...");

    // at this point we have bigL and bigW all set

    LOG("Computing eigenvalues and eigenvectors by power iteration...");

    // double prev_lambda = 1.0;
    // double lambda = 0.0;
    uint iter_num = 0;
    uint num_eigens = 10;
    // arma::vec prev_x = arma::randn<arma::vec>(num_pixels+len_spectrum);
    arma::mat prev_x = arma::randn<arma::mat>(num_pixels+len_spectrum, num_eigens);
    arma::mat x = bigL * prev_x;
    arma::mat q, r_mat;
    arma::vec lambda = arma::randn<arma::vec>(num_eigens);
    arma::vec prev_lambda = arma::randn<arma::vec>(num_eigens);
    while (arma::norm(lambda - prev_lambda) > 0.001) {
        ++iter_num;
        prev_lambda = lambda;
        prev_x = x;
        x = bigL * x;
        for (uint i=0; i<num_eigens; ++i) {
            lambda(i) = x(0,i) / prev_x(0,i);
        }
        arma::qr_econ(q, r_mat, x);
        x = q;
        LOG("\t\t" << iter_num << "\tdifflambda=" << norm(lambda - prev_lambda) << "lambdas:\n" << lambda);
        // LOG("\t\t\t" << lambda << "\t" << arma::norm(bigL*x - lambda*x) << "\t" << arma::norm(prev_x - x));
    }
    for (uint i=0; i < num_eigens; ++i) {
        LOG("\teigenvalue " << i << " = " << lambda(i) << " with error " << arma::norm(bigL*x.col(i) - lambda(i)*x.col(i)));
    }
    
    // x.print();
    // LOG("\tcheck against ");
    // arma::vec test_x = (1/lambda)*(bigL * x);
    // test_x.print();

    // REDSVD::RedSVD rsvd(bigL, 10);


    // LOG("Computing SVD...");

    // REDSVD::RedSVD rsvd(bigL, 10);



    // arma::cx_vec eigval;
    // arma::cx_mat eigvec;

    // arma::vec s = svd(bigL);

    // eig_gen(eigval, eigvec, bigL);

    // uint mat_size = num_pixels + len_spectrum;
    // gsl_eigen_gen_workspace * gsl_workspace =  gsl_eigen_gen_alloc (mat_size);
    // gsl_vector_complex * vec_alpha = gsl_vector_complex_alloc(mat_size);
    // gsl_vector * vec_beta = gsl_vector_alloc(mat_size);

    // gsl_eigen_gen (bigL, bigW, vec_alpha, vec_beta, gsl_workspace);

    // printf("solution:\n");
    // {
    //   gsl_vector_view a_real = gsl_vector_complex_real (vec_alpha);
    //   gsl_vector_view a_imag = gsl_vector_complex_imag (vec_alpha);
    //   int k;
    //   for (k = 0; k < 4; k++) { 
    //     printf("% .18e % .18e % .18e\n", 
    //            gsl_vector_get(&a_real.vector,k), 
    //            gsl_vector_get(&a_imag.vector,k),
    //            gsl_vector_get(vec_beta,k));
    //   }
    // }

    // gsl_eigen_gen_free (gsl_workspace);


    // gsl_matrix_view g_bigW = gsl_matrix_view_array (bigW, num_pixels + len_spectrum, num_pixels + len_spectrum);

    delete maxima;
    delete specsdiag;
    delete pixdiag;
    delete on_diag;
    delete little;
    delete_2d<double>(len_spectrum, spectra);
    delete_2d<double>(num_pixels, W);

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

