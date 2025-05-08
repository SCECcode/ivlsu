/**
 * @file test.c
 * @brief Bootstraps the test framework for the IMPERIAL/IVLSU library.
 * @author - SCEC
 * @version 1.0
 *
 * Tests the IMPERIAL/IVLSU library by loading it and executing the code as
 * UCVM would do it.
 *
 */

#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include "ivlsu.h"

/**
 * Initializes and runs the test program. Tests link against the
 * static version of the library to prevent any dynamic loading
 * issues.
 *
 * @param argc The number of arguments.
 * @param argv The argument strings.
 * @return A zero value indicating success.
 */
int main(int argc, const char* argv[]) {

	// Declare the structures.
	ivlsu_point_t pt;
	ivlsu_properties_t ret;

	// Initialize the model. 
        // try to use Use UCVM_INSTALL_PATH
        char *envstr=getenv("UCVM_INSTALL_PATH");
        if(envstr != NULL) {
	   assert(ivlsu_init(envstr, "ivlsu") == 0);
           } else {
	     assert(ivlsu_init("..", "ivlsu") == 0);
        }

	printf("Loaded the model successfully.\n");

	// Query a point.
	pt.longitude = -116.0516;
	pt.latitude = 32.6862;
	pt.depth = 2000;

	ivlsu_query(&pt, &ret, 1);

        printf("vs : %lf\n",ret.vs);
        printf("vp : %lf\n",ret.vp);
        printf("rho: %lf\n",ret.rho);

	assert(ret.vs = 2906.943265);
	assert(ret.vp = 4838.500000);
	assert(ret.rho = 2510.425129);

	printf("Query was successful.\n");

	// Close the model.
	assert(ivlsu_finalize() == 0);

	printf("Model closed successfully.\n");

	printf("\nALL IMPERIAL TESTS PASSED\n");

	return 0;
}
