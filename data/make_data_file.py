#!/usr/bin/env python


##
#  Builds the data files in the expected format from IV33.dat.txt
#
# from >> easting(Km) northing(Km) depth(Km) vp(km/s)
#
# to >>  /** P-wave velocity in km/s per second */
#        double vp;
# depth is in increment of 1000m,
#
#The columns of the file are: x (km) y (km) z (km) vp (km/s), where 
#x, y, and z are utm coordinates in km and the columns increase in x.y.z order. 
#Grid spacing is 1 km in each spatial dimension. 

import getopt
import sys
import subprocess
import struct
import array
import pdb

##import osr
##def transform_utm_to_wgs84(easting, northing, zone):
##    utm_coordinate_system = osr.SpatialReference()
##          # Set geographic coordinate system to handle lat/lon
##    utm_coordinate_system.SetWellKnownGeogCS("WGS84") 
##
##    is_northern = northing > 0    
##    utm_coordinate_system.SetUTM(zone, is_northern)
##          # Clone ONLY the geographic coordinate system 
##    wgs84_coordinate_system = utm_coordinate_system.CloneGeogCS() 
##          # create transform component
##          # (<from>, <to>)
##    utm_to_wgs84_transform = osr.CoordinateTransformation(utm_coordinate_system, wgs84_coordinate_system) 
##          # returns lon, lat, altitude
##    return utm_to_wgs84_transform.TransformPoint(easting, northing, 0) 


## at hypocenter  LSU/IV33.dat.txt

model = "LSU"

dimension_x = 66
dimension_y = 86 
dimension_z = 9

lon_origin = -116.051578
lat_origin = 32.596922

lon_upper = -115.344866
lat_upper = 33.356203

delta_lon = (lon_upper - lon_origin )/(dimension_x-1)
delta_lat = (lat_upper - lat_origin)/(dimension_y-1)

def usage():
    print("\n./make_data_files.py -u [uid]\n\n")
    print("-u - username to use to do the dataset retrieval.\n")
    sys.exit(0)

def main():

    # Set our variable defaults.
    username = ""
    path = "/var/www/html/research/ucvmc/" + model 

    try:
        opts, args = getopt.getopt(sys.argv[1:], "u:", ["user="])
    except getopt.GetoptError as err:
        print(str(err))
        usage()
        sys.exit(1)

    for o, a in opts:
        if o in ("-u", "--user"):
            username = str(a) + "@"

##
##    subprocess.check_call(["scp", username +
##                           "hypocenter.usc.edu:" + path + "/IV33.dat.txt",
##                           "."])
##

    # Now we need to go through the data files and put them in the correct
    # format for LSU_IV. More specifically, we need a Vp.dat

    f = open("./IV33.dat.txt")

    f_vp = open("./iv/vp.dat", "wb")
    f_easting = open("./iv/easting.dat", "wb")
    f_northing = open("./iv/northing.dat", "wb")

    vp_arr = array.array('f', (-1.0,) * (dimension_x * dimension_y * dimension_z))
    easting_arr = array.array('i', (-1,) * (dimension_x * dimension_y * dimension_z))
    northing_arr = array.array('i', (-1,) * (dimension_x * dimension_y * dimension_z))

    print "dimension is", (dimension_x * dimension_y * dimension_z)

    nan_cnt = 0
    total_cnt =0;
    x_pos=0;
    y_pos=0;
    z_pos=0;
    for line in f:
        arr = line.split()

        easting_v = int(arr[0])
        northing_v = int(arr[1])
        depth_v = int(arr[2])
        vp_v = float(arr[3])
        total_cnt = total_cnt + 1

        vp = vp_v * 1000.0;
        loc =z_pos * (dimension_y * dimension_x) + (y_pos * dimension_x) + x_pos
        vp_arr[loc] = vp
        easting_arr[loc] = easting_v
        northing_arr[loc] = northing_v

        print total_cnt, "loc",loc," ", x_pos," ",y_pos," ",z_pos," >> ",easting_v," ",northing_v," ",depth_v,":",vp 

      
        x_pos = x_pos + 1
        if(x_pos == dimension_x) :
          x_pos = 0;
          y_pos = y_pos+1
          if(y_pos == dimension_y) :
            y_pos=0;
            z_pos = z_pos+1
            if(z_pos == dimension_z) :
              print "All DONE"

    vp_arr.tofile(f_vp)
    easting_arr.tofile(f_easting)
    northing_arr.tofile(f_northing)

    f.close()
    f_vp.close()
    f_easting.close()
    f_northing.close()

    print("Done! total", total_cnt)

if __name__ == "__main__":
    main()

