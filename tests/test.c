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
	assert(ivlsu_init("../", "ivlsu") == 0);

	printf("Loaded the model successfully.\n");

	// Query a point.
	pt.longitude = -118;
	pt.latitude = 34;
	pt.depth = 0;

	ivlsu_query(&pt, &ret, 1);

	assert(ret.vs > 0);
	assert(ret.vp > 0);
	assert(ret.rho > 0);

	printf("Query was successful.\n");

	// Close the model.
	assert(ivlsu_finalize() == 0);

	printf("Model closed successfully.\n");

	printf("\nALL IMPERIAL TESTS PASSED");

	return 0;
}
