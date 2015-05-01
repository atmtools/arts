/*!
  \file   test_oem.cc
  \author Simon Pfreundschuh <simon@thinks>
  \date   Sat Apr 18 19:50:30 2015

  \brief  Test for the OEM functions.
*/


#include <cmath>
#include "engine.h"
#include <fstream>
#include "lin_alg.h"
#include <iostream>
#include <iomanip>
#include "matrix.h"
#include "oem.h"
#include <stdlib.h>
#include <string>
#include <time.h>
#include "unistd.h"

using std::abs;
using std::cout;
using std::endl;
using std::min;
using std::ofstream;
using std::setw;
using std::string;

string source_dir = SOURCEDIR;
string atmlab_dir = ATMLABDIR;

//! Fill matrix with random values.
/*!

Fills the given matrix with integer values in the range [0, range] or
[-range, range], if positive is true.

  \param[out] A The matrix to be filled.
  \param[in] range The range of the values to fill the matrix with.
  \param positive If true the matrix is filled with values from the interval
                  [0,range], otherwise the values are taken from the interval
		  [-range, range].
*/
void random_fill_matrix( MatrixView A,
			 Index range = RAND_MAX,
			 bool positive = true )
{
    Index m = A.nrows();
    Index n = A.ncols();

    for (Index i=0; i<m; i++)
    {
	for (Index j=0; j<n; j++)
	{
	    if (positive)
	    {
		A(i,j) = (Numeric) (rand() % range);
	    } else {
		A(i,j) = (Numeric) (rand() % (2 * range) - range);
	    }
	}
    }
}

//! Fill vector with random values.
/*!

  Fills the given vector with random integer values from the range [0, range], or,
  if positive is set to true, from the range [-range, range].

  \param[out] v The vector to be filled.
  \param[in] range The range from which the values are taken.
  \param[in] positive If true, the values are taken from the interval [0, range],
                      otherwise from the range [-range, range].
*/
void random_fill_vector( VectorView v,
			 Index range = RAND_MAX,
			 bool positive = true )
{
    Index n = v.nelem();
    for (Index i = 0; i < n; i++)
    {
	if ( positive )
	{
	    v[i] = (Numeric) (rand() % range);
	} else {
	    v[i] = (Numeric) (rand() % (2 * range) - range);
	}
    }
}

//! Write matrix to text file.
/*!

  Write the given matrix in plain text to the file filename in the
  current directory.
  \param[in] A The matrix to write to the file.
  \param[in] filename The name of the file to write to.
*/
void write_matrix( ConstMatrixView A,
		   const char* filename )
{

    Index m = A.nrows();
    Index n = A.ncols();

    ofstream ofs( filename, ofstream::out);

    for (Index i = 0; i < m; i++)
    {
	for (Index j = 0; j < (n - 1); j++)
	{
	    ofs << A(i,j) << " ";
	}
	ofs << A( i, n - 1 );
	ofs << endl;
    }
    ofs.close();
}

//! Write vector to text file.
/*!
  Write the given vector to the file filename in the current directory.
  \param[in] v The vector to write.
  \param[in] filename The name of the file to write to.
*/void write_vector( ConstVectorView v,
		   const char* filename )
{
    Index n = v.nelem();

    ofstream ofs( filename, ofstream::out);

    for ( Index i=0; i<n; i++)
    {
	ofs << v[i] << endl;
    }
    ofs.close();
}

//! Fill matrix with the given value.
/*!

  \param[in,out] A The matrix to be filled.
  \param[in] value The value to fill the matrix with.
*/
void fill_matrix( MatrixView A,
		  Numeric value )
{

    Index m = A.nrows();
    Index n = A.ncols();

    for (Index i=0; i<m; i++)
    {
	for (Index j=0; j<n; j++)
	{
	    A( i, j ) = value;
	}
    }
}

//! Fill vector with the given value.
/*!

  \param[in,out] v The vector to be filled.
  \param[in] value The value to fill the vector with.
*/
void fill_vector( VectorView v,
		  Numeric value )
{
    Index n = v.nelem();

    for ( Index i=0; i<n; i++)
    {
	v[i] = value;
    }
}

//! Generate test data for linear OEM retrieval.
/*!
  Fills the given matrices and vectors needed for the linear OEM retrieval
  functions with random values.

  \param y The measurement vector. Filled with random integer values in the range
           [0, 10].
  \param xa The a priori vector. Filled with random integer values in the range
           [0, 10].
  \param K The linear forward model. Filled with random integer values in the
           range [-10, 10].
  \param Se The covariance matrix for observation uncertainties. Filled with
            random values in the range [-1, 1] (Scaled down from the range
	    [-10, 10]).
  \param Sa The covariance matrix for a priori uncertainties. Filled with random
            values in the range [-1, 1] (Scaled down from the range [-10, 10]).
*/
void generate_test_data( VectorView y,
			 VectorView xa,
			 MatrixView K,
			 MatrixView Se,
			 MatrixView Sa)
{
    random_fill_vector( y, 10 );
    random_fill_vector( xa, 10 );
    random_fill_matrix( K, 10, false );

    random_fill_matrix( Se, 10, false);
    Se *= 0.1;

    random_fill_matrix( Sa, 10, false);
    Sa *= 0.1;
}

//! Maximum element-wise error of two vectors.
/*!
  If relative is true, the maximum element-wise error is computed.
  Otherwise the absolute error is computed.

  \param[in] v1 The first vector.
  \param[in] v2 The reference vector used to normalize the relative error.
  \param relative If true the relative error is computed, otherwise the absolute
                  error is computed.

  \return The maximum relative or absolute element-wise error.
*/
Numeric max_error_vector( ConstVectorView v1,
			  ConstVectorView v2,
			  bool relative = true )
{
    Index n = min( v1.nelem(), v2.nelem() );
    Numeric max = 0.0, err = 0.0;
    for ( Index i = 0; i < n; i++ )
    {
	err = 0.0;

	if (relative)
	{

	    if (v2[i] != 0.0)
	    {
		err = abs( v2[i] - v1[i] ) / abs( v2[i] );
	    }

	} else {

	    err = abs( v2[i] - v1[i] );

	}

	if (err > max)
	    max = err;
    }
    return err;
}

//! Run test script in matlab.
/*!
 Runs the test script given by filename in matlab. Reads out the
 the execution time from the variable t in the current workspace and
 returns is as return value.

  \param[in] eng Pointer to the running matlab engine.
  \param[in] filename Name of the file to be run.

  \return The time required for the execution.
*/
Index run_test_matlab( Engine* eng,
                       string filename )
{

    mxArray* t;
    Index time;

    // Run test.
    string cmd = "run('" + filename + "');";
    engEvalString( eng, cmd.c_str() );

    // Get execution time from matlab.
    t = engGetVariable( eng, "t" );
    time = (Index) ((Numeric*) mxGetData( t ))[0];
    return time;

}

//! Run linear OEM in matlab.
/*!
  Runs the oem function from the atmlab package. Runs the external matlab
  script "test_oem.m". The test data is read from the text files "y.txt",
  "x_a.txt", "K.txt", "Se.txt", "Sa.txt" from the current directory. Writes
  the results into the provided vector. Returns the execution time measured
  in matlab with the 'cputime' function.

  \param[out] x Vector to write the results of the retrieval to.
  \param[in] eng Pointer to the matlab engine that manages the running Matlab
                 session.

  \return Execution time of oem in matlab in ms.
*/
Index run_oem_matlab( VectorView x,
		      Engine* eng )
{
    Index n = x.nelem(), time;
    mxArray *x_m, *t;

    // Run test.
    engEvalString( eng, "run('test_oem')");

    // Read out results.
    x_m = engGetVariable( eng, "x" );

    for ( Index i = 0; i < n; i++ )
    {
	x[i] = ((Numeric*) mxGetData( x_m ))[i];
    }

    // Get execution time from matlab.
    t = engGetVariable( eng, "t" );
    time = (Index) ((Numeric*) mxGetData( t ))[0];
    return time;
}


//! Setup the test environment.
/*!
  Changes to the test directory and initializes the matlab engine. Initialized
  the atmlab package.

  \param[in,out] eng Pointer variable that will contain the pointer to the
                     initialized Matlab engine.
*/
void setup_test_environment( Engine * &eng )
{
    // swith to test folder
    string cmd;
    cmd = source_dir + "/test_oem_files";
    int out = chdir( cmd.c_str() );
    (void) out;

    // Start MATLAB and try to initialize atmlab package.
    string atmlab_init = "run('" + atmlab_dir + "/atmlab/atmlab_init.m');";

    eng = engOpen(NULL);

    engEvalString( eng, atmlab_init.c_str() );
    cmd = "cd('" + source_dir + "/test_oem_files');";
    engEvalString( eng, cmd.c_str() );
}

//! Plot benchmark results
/*!

  Run matlab script that generates a plot of the benchmark results.

  \param[in] eng Pointer to the running matlab engine.
  \param[in] filename Filename of the file containing the results. Should be
                      'times_mult.txt' or 'times_linear.txt'.
*/
void run_plot_script( Engine *eng,
		      string filename,
                      string title )
{

    string cmd = "filename = '" + filename + "'";
    engEvalString( eng, cmd.c_str() );
    cmd = "plot_title = '" + title + "'";
    engEvalString( eng, cmd.c_str() );
    engEvalString( eng, "run('make_plot.m');" );

}

//! Tidy up test environment
/*!
  Deletes temporary test files and closes the Matlab session.
  \param[in] eng Pointer to the running Matlab engine.
*/
void tidy_up_test_environment( Engine *eng)
{
    int out = system( "rm *_t.txt" );
    (void) out;

    engEvalString( eng, "close()" );
}

//! Matrix inversion benchmark.
/*!

  Inverts randomly generated matrix matrices in matlab and in arts and compares
  the performance. Performs ntests numbers of test with matrices linearly
  increasing in size starting at n0 and ending an n1. Writes the result to
  standard out and to the file "times_inv.txt" in the current directory. Also
  generates a plot of the data using matlab stored in "times_inv.png".

  \param[in] n0 Size of the smallest matrix in the benchmark.
  \param[in] n1 Size of the largest matrix in the benchmark.
  \param[in] ntests Number of tests to be performed. In each step the size
                    of the matrix is linearly increased from n0 to n1.
*/
void benchmark_inv( Engine* eng,
		    Index n0,
		    Index n1,
		    Index ntests )
{
    Index step = (n1 - n0) / (ntests - 1);
    Index n = n0;

    ofstream ofs( "times_inv.txt", ofstream::out );
    ofs << "#" << setw(4) << "n" << setw(10) << "C++";
    ofs << setw(10) << "Matlab" << endl;

    cout << endl << "N TIMES N MATRIX INVERSION" << endl << endl;
    cout << setw(5) << "n" << setw(10) << "C++" << setw(10);
    cout << "Matlab" << endl;

    for ( Index i = 0; i < ntests; i++ )
    {
	Matrix A(n,n), B(n,n);

	random_fill_matrix( A, 100 );
	write_matrix( A, "A_t.txt");

	Index t, t1, t2, t_m;

	t1 = clock();
	inv( B, A );
	t2 = clock();
	t = (t2 - t1) * 1000 / CLOCKS_PER_SEC;
	t_m = run_test_matlab( eng, "test_inv.m" );

	ofs << setw(5) << n << setw(10) << t << setw(10) << t_m << endl;
	cout << setw(5) << n << setw(10) << t << setw(10) << t_m << endl;

	n += step;
    }
    cout << endl << endl;

    // Tidy up
    ofs.close();
    run_plot_script( eng, "times_inv.txt", "Matrix Inversion" );

}

//! Matrix multiplication benchmark.
/*!

  Multiplies two identical matrices in matlab and in arts and compares
  the performance. Performs ntests numbers of test with matrices linearly
  increasing in size starting at n0 and ending an n1. Writes the result to
  standard out and to the file "times_mult.txt" in the current directory. Also
  generates a plot of the data using matlab stored in "times_mult.png".

  \param[in] n0 Size of the smallest matrix in the benchmark.
  \param[in] n1 Size of the largest matrix in the benchmark.
  \param[in] ntests Number of tests to be performed. In each step the size
                    of the matrix is linearly increased from n0 to n1.
*/
void benchmark_mult( Engine* eng,
		     Index n0,
		     Index n1,
		     Index ntests )
{
    Index step = (n1 - n0) / (ntests - 1);
    Index n = n0;

    ofstream ofs( "times_mult.txt", ofstream::out );
    ofs << "#" << setw(4) << "n" << setw(10) << "Matlab";
    ofs << setw(10) << "C++" << endl;

    cout << endl << "N TIMES N MATRIX MULTIPLICATION" << endl << endl;
    cout << setw(5) << "n" << setw(10) << "C++" << setw(10);
    cout << "Matlab" << endl;

    for ( Index i = 0; i < ntests; i++ )
    {
	Matrix A(n,n), B(n,n);

	random_fill_matrix( A, 100 );
	write_matrix( A, "A_t.txt");

	Index t, t1, t2, t_m;

	t1 = clock();
	mult( B, A, A );
	t2 = clock();
	t = (t2 - t1) * 1000 / CLOCKS_PER_SEC;
	t_m = run_test_matlab( eng, "test_mult.m" );

	ofs << setw(5) << n << setw(10) << t << setw(10) << t_m << endl;
	cout << setw(5) << n << setw(10) << t << setw(10) << t_m << endl;

	n += step;
    }
    cout << endl << endl;

    // Tidy up
    ofs.close();
    run_plot_script( eng, "times_mult.txt", "Matrix Multiplication" );

}

//! Benchmark linear oem.
/*!
  Run linear oem test and benchmark. Runs ntests numbers of tests with the size
  of K growing linearly from n0 to n1. Prints for each test the maximum relative
  error in the result compared to the matlab implementation and the cpu time in
  milliseconds needed for the execution. The timing results are also written to
  to a file "times_linear.txt" in the current directory.

  \param n0 Starting size for the K matrix.
  \param n1 Final size for the K matrix.
  \param ntests Number of tests to be performed.
*/
void benchmark_linear_oem( Engine* eng,
			   Index n0,
			   Index n1,
			   Index ntests )
{
    Index step = (n1 - n0) / (ntests - 1);
    Index n = n0;

    ofstream ofs( "times_linear.txt", ofstream::out );

    ofs << "#" << setw(4) << "n" << setw(10) << "Matlab";
    ofs << setw(10) << "C++" << endl;

    cout << endl << "LINEAR OEM" << endl << endl;
    cout << setw(5) << "n" << setw(10) << "C++" << setw(10);
    cout << "Matlab" << setw(20) << "Max. Rel. Error" << endl;


    // Run tests.
    for ( Index i = 0; i < ntests; i++ )
    {
	Vector x(n), x_m(n), y(n), xa(n);
	Matrix K(n,n), Se(n,n), Sa(n,n);

	generate_test_data( y, xa, K, Se, Sa );

	write_vector( xa, "xa_t.txt" );
	write_vector( y, "y_t.txt" );
	write_matrix( K, "K_t.txt" );
	write_matrix( Se, "Se_t.txt" );
	write_matrix( Sa, "Sa_t.txt" );

	Index t, t1, t2, t_m;

	t1 = clock();
	linear_oem( x, y, xa, K, Se, Sa, false );
	t2 = clock();
	t = (t2 - t1) * 1000 / CLOCKS_PER_SEC;
	t_m = run_oem_matlab( x_m, eng );

	ofs << setw(5) << n << setw(10) << t << setw(10) << t_m << endl;
	cout << setw(5) << n << setw(10) << t << setw(10) << t_m;
	cout << setw(20) << max_error_vector( x, x_m ) << endl;

	n += step;
    }
    cout << endl << endl;

    // Tidy up
    ofs.close();
    run_plot_script( eng, "times_linear.txt", "Linear OEM" );

}

//! Test linear_oem.
/*!
  Tests the linear_oem function using randomized input data. Performs
  ntest numbers of tests with a m times n K-matrix. For each test, the
  maximum relative error is printed to standard out.

  \param[in] m Size of the measurement space.
  \param[in] n Size of the state space.
  \param[in] ntests Number of tests to be performed.
*/
void test_linear_oem( Engine* eng,
		      Index m,
		      Index n,
		      Index ntests )
{
    Vector x(n), x_m(n), y(m), xa(n);
    Matrix K(m,n), Se(m,m), Sa(n,n);

    // Run tests.
    for ( Index i = 0; i < ntests; i++ )
    {
	generate_test_data( y, xa, K, Se, Sa );

	write_vector( xa, "xa_t.txt" );
	write_vector( y, "y_t.txt" );
	write_matrix( K, "K_t.txt" );
	write_matrix( Se, "Se_t.txt" );
	write_matrix( Sa, "Sa_t.txt" );

	linear_oem( x, y, xa, K, Se, Sa, true );
	write_vector( x, "x_t.txt" );

	run_oem_matlab( x_m, eng );
	cout << max_error_vector( x, x_m ) << endl;
    }
}

int main()
{

    // Set up the test environment.
    Engine * eng;
    setup_test_environment( eng );

    // Run tests and benchmarks.
    benchmark_inv( eng, 100, 1000, 10);
    benchmark_mult( eng, 100, 1000, 10);
    benchmark_linear_oem( eng, 100, 1000, 10);

    // Tidy up test environment.
    tidy_up_test_environment( eng );

}
