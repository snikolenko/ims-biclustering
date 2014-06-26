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

int help(char *argv[]) {
    cout << "Usage: " << argv[0] << " [--mat] --input=input_file" << endl;
    return 0;
}

int main(int argc, char *argv[]) {
    const struct option longopts[] = {
        {"help",            no_argument,        0, 'h'},
        {"input",           required_argument,  0, 'i'},
        {"mat",             no_argument,        0, 'm'},
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
            case 'm':   matlab_input = true;            break;
        }
    }

    mat_t *matfp;
    matvar_t *matvar;

    if (matlab_input) {
        LOG("Reading input from " << input_filename << " as a Matlab file...");
        matfp = Mat_Open(input_filename.c_str(), MAT_ACC_RDONLY);

        if ( NULL == matfp ) {
            fprintf(stderr,"Error creating MAT file \"matfile5.mat\"!\n");
            return EXIT_FAILURE;
        }

        LOG("Variables in " << input_filename << ":");
        while ( (matvar = Mat_VarReadNextInfo(matfp)) != NULL ) {
            LOG("\t" << matvar->name);
            Mat_VarFree(matvar);
            matvar = NULL;
        }
    }

    double data[] = { 1.0  , 1/2.0, 1/3.0, 1/4.0,
                    1/2.0, 1/3.0, 1/4.0, 1/5.0,
                    1/3.0, 1/4.0, 1/5.0, 1/6.0,
                    1/4.0, 1/5.0, 1/6.0, 1/7.0 };

    gsl_matrix_view m = gsl_matrix_view_array (data, 4, 4);

    gsl_vector *eval = gsl_vector_alloc (4);
    gsl_matrix *evec = gsl_matrix_alloc (4, 4);

    gsl_eigen_symmv_workspace * w =  gsl_eigen_symmv_alloc (4);
    gsl_eigen_symmv (&m.matrix, eval, evec, w);
    gsl_eigen_symmv_free (w);
    gsl_eigen_symmv_sort (eval, evec, GSL_EIGEN_SORT_ABS_ASC);

    LOG("Sample output from gsl_eigen:")  
    for (uint i = 0; i < 4; ++i) {
        double eval_i = gsl_vector_get (eval, i);
        gsl_vector_view evec_i = gsl_matrix_column (evec, i);

        LOG ("\teigenvalue = " << eval_i);
        LOG ("\teigenvector = ");
        gsl_vector_fprintf (stdout, &evec_i.vector, "%g");
    }

    LOG("All done.");
    if (matlab_input) {
        Mat_Close(matfp);
    }
    return 0;
}
