#!/usr/bin/env python

##
#  Retrieve the easting/norhting from the original data txt file and
#  and retrieve the material properties from binary vp.dat file 
#
## import utm
## (rlat_v, rlon_v) = utm.to_latlon(easting_v*1000, northing_v*1000, zone, 'S')
## lat_v=float('%.6f'%rlat_v)
## lon_v=float('%.6f'%rlon_v)

import getopt
import sys
import subprocess
import struct
import numpy as np

dimension_x = 66
dimension_y = 86 
dimension_z = 9

lon_origin = -116.051578
lat_origin = 32.596922
lon_upper = -115.344866
lat_upper = 33.356203

easting_origin = 589
northing_origin = 3607
easting_upper = 654
northing_upper = 3692

delta_easting = (easting_upper - easting_origin )/(dimension_x-1)
delta_northing = (northing_upper - northing_origin)/(dimension_y-1)

def usage():
    print("\n./query_data_files.py\n\n")
    sys.exit(0)

def main():
    total_cnt=0
    f = open("./IV33.dat.txt")
    f_vp = open("./iv/vp.dat")
    vp_arr = np.fromfile(f_vp, dtype=np.float32)

    x_pos=0;
    y_pos=0;
    z_pos=0;
    for line in f:
        arr = line.split()

        easting_v = int(arr[0])
        northing_v = int(arr[1])
        depth_v = int(arr[2])
        rvp_v = float(arr[3]) * 1000.0
        vp_v = float("%0.1f" % rvp_v)

        x_pos = int(round((easting_v - easting_origin) / delta_easting))
        y_pos = int(round((northing_v - northing_origin) / delta_northing))
        z_pos = int(depth_v)

        loc =z_pos * (dimension_y * dimension_x) + (y_pos * dimension_x) + x_pos
        fvp= float("%0.1f" % vp_arr[loc].item())
        
        if(vp_v != fvp):
           print total_cnt,"BAD",vp_v,"fvp",fvp
           
        total_cnt = total_cnt + 1

    f_vp.close()
    print "DONE with total_cnt ",total_cnt

## dump all
def dump():
    f_vp = open("./iv/vp.dat")
    f_vp.close()
    count=0
    x_pos =0
    y_pos =0
    z_pos =0

    while(1) :
      offset=z_pos * (dimension_y * dimension_x) + (y_pos * dimension_x) + x_pos
      vp=vp_arr[offset];
      print count,"xyz:", x_pos,y_pos,z_pos,">> vp", vp

      count = count + 1
      x_pos = x_pos + 1
      if(x_pos == dimension_x) :
        x_pos = 0;
        y_pos = y_pos+1
        if(y_pos == dimension_y) :
          y_pos=0;
          z_pos = z_pos+1
          if(z_pos == dimension_z) :
            print "Done! count ",count
            exit(0)


if __name__ == "__main__":
    main()


