/**
 * @file imperial.h
 * @brief Main header file for LSU Coachella Valley library.
 * @author - SCEC 
 * @version 1.0
 *
 * Delivers LSU Coachella Valley Velocity Model
 *
 */

// Includes
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <math.h>

#include "proj_api.h"

// Constants
#ifndef M_PI
	/** Defines pi */
	#define M_PI 3.14159265358979323846
#endif

/** Defines a return value of success */
#define SUCCESS 0
/** Defines a return value of failure */
#define FAIL 1
/** Defines a return value of NA from model */
#define NA -1 

// Structures
/** Defines a point (latitude, longitude, and depth) in WGS84 format */
typedef struct imperial_point_t {
	/** Longitude member of the point */
	double longitude;
	/** Latitude member of the point */
	double latitude;
	/** Depth member of the point */
	double depth;
} imperial_point_t;

/** Defines the material properties this model will retrieve. */
typedef struct imperial_properties_t {
	/** P-wave velocity in meters per second */
	double vp;
	/** S-wave velocity in meters per second */
	double vs;
	/** Density in g/m^3 */
	double rho;
} imperial_properties_t;

/** The IMPERIAL configuration structure. */
typedef struct imperial_configuration_t {
	/** The zone of UTM projection */
	int utm_zone;
	/** The model directory */
	char model_dir[128];
	/** Number of x points */
	int nx;
	/** Number of y points */
	int ny;
	/** Number of z points */
	int nz;
	/** Depth in meters */
	double depth;
	/** Top left corner easting */
	double top_left_corner_e;
	/** Top left corner northing */
	double top_left_corner_n;
	/** Top right corner easting */
	double top_right_corner_e;
	/** Top right corner northing */
	double top_right_corner_n;
	/** Bottom left corner easting */
	double bottom_left_corner_e;
	/** Bottom left corner northing */
	double bottom_left_corner_n;
	/** Bottom right corner easting */
	double bottom_right_corner_e;
	/** Bottom right corner northing */
	double bottom_right_corner_n;
	/** Z interval for the data */
	double depth_interval;
	double p5;
        /** Bilinear or Trilinear Interpolation on or off (1 or 0) */
        int interpolation;

} imperial_configuration_t;

/** The model structure which points to available portions of the model. */
typedef struct imperial_model_t {
	/** A pointer to the Vp data either in memory or disk. Null if does not exist. */
	void *vp;
	/** Vp status: 0 = not found, 1 = found and not in memory, 2 = found and in memory */
	int vp_status;
} imperial_model_t;

// Constants
/** The version of the model. */
const char *imperial_version_string = "IMPERIAL";

// Variables
/** Set to 1 when the model is ready for query. */
int imperial_is_initialized = 0;

/** Location of the binary data files. */
char imperial_data_directory[128];

/** Configuration parameters. */
imperial_configuration_t *imperial_configuration;
/** Holds pointers to the velocity model data OR indicates it can be read from file. */
imperial_model_t *imperial_velocity_model;

/** Proj.4 latitude longitude, WGS84 projection holder. */
projPJ imperial_latlon;
/** Proj.4 UTM projection holder. */
projPJ imperial_utm;

/** The cosine of the rotation angle used to rotate the box and point around the bottom-left corner. */
double imperial_cos_rotation_angle = 0;
/** The sine of the rotation angle used to rotate the box and point around the bottom-left corner. */
double imperial_sin_rotation_angle = 0;

/** The height of this model's region, in meters. */
double imperial_total_height_m = 0;
/** The width of this model's region, in meters. */
double imperial_total_width_m = 0;

// UCVM API Required Functions

#ifdef DYNAMIC_LIBRARY

/** Initializes the model */
int model_init(const char *dir, const char *label);
/** Cleans up the model (frees memory, etc.) */
int model_finalize();
/** Returns version information */
int model_version(char *ver, int len);
/** Queries the model */
int model_query(imperial_point_t *points, imperial_properties_t *data, int numpts);

#endif

// IMPERIAL Related Functions

/** Initializes the model */
int imperial_init(const char *dir, const char *label);
/** Cleans up the model (frees memory, etc.) */
int imperial_finalize();
/** Returns version information */
int imperial_version(char *ver, int len);
/** Queries the model */
int imperial_query(imperial_point_t *points, imperial_properties_t *data, int numpts);

// Non-UCVM Helper Functions
/** Reads the configuration file. */
int imperial_read_configuration(char *file, imperial_configuration_t *config);
void print_error(char *err);
/** Retrieves the value at a specified grid point in the model. */
void imperial_read_properties(int x, int y, int z, imperial_properties_t *data);
/** Attempts to malloc the model size in memory and read it in. */
int imperial_try_reading_model(imperial_model_t *model);
/** Calculates density from Vs. */
double imperial_calculate_density(double vs);

// Interpolation Functions
/** Linearly interpolates two imperial_properties_t structures */
void imperial_linear_interpolation(double percent, imperial_properties_t *x0, imperial_properties_t *x1, imperial_properties_t *ret_properties);
/** Bilinearly interpolates the properties. */
void imperial_bilinear_interpolation(double x_percent, double y_percent, imperial_properties_t *four_points, imperial_properties_t *ret_properties);
/** Trilinearly interpolates the properties. */
void imperial_trilinear_interpolation(double x_percent, double y_percent, double z_percent, imperial_properties_t *eight_points,
							 imperial_properties_t *ret_properties);
