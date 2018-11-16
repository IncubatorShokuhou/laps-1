#!/bin/sh 

# Annotate a simulated weather image
# This script should run in the same directory as the label files
# Variables $ILOC, $MODE_ALLSKY, $POINT, and $YDISP are inherited from environment

  echo "start annotate_allsky.sh for $ILOC"
  pwd

  echo "MODE_ALLSKY = $MODE_ALLSKY"
  echo "POINT = $POINT"

  ALT_SCALE=`head -3 label2.$ILOC | tail -1 | cut -c17-23`
  AZI_SCALE=`head -3 label2.$ILOC | tail -1 | cut -c24-30`
  LATLON=`head -1 label2.$ILOC`
  LATLON=`echo $LATLON | sed 's/^[ \t]*//'` # remove leading spaces
  HT=`head -5 label2.$ILOC | tail -1`
  CORRELATION=`head -7 label2.$ILOC | tail -1 | awk '{print $4}'`

  echo "AZI_SCALE = $AZI_SCALE"
  echo "CORRELATION = $CORRELATION"

# Annotate Model
  if test "$MODE_ALLSKY" = "polar" || test "$MODE_ALLSKY" = "both"; then
    IMGGEOM=`identify allsky_polar_$ILOC.png | awk '{print tolower($3)}'`
    echo "polar IMGGEOM = $IMGGEOM"

    if test "$IMGGEOM" = "511x511"; then
      POINTP=16
    elif test "$IMGGEOM" = "1023x1023"; then
      POINTP=22
    elif test "$IMGGEOM" = "1535x1535"; then
      POINTP=22
    elif test "$IMGGEOM" = "2047x2047"; then
      POINTP=22
    elif test "$IMGGEOM" = "2559x2559"; then
      POINTP=22
    else
      POINTP=22
    fi

    convert -fill white -annotate +5+20  "Simulated" -pointsize $POINTP allsky_polar_$ILOC.png allsky_polar_$ILOC.png
  fi

  if test "$MODE_ALLSKY" = "cyl" || test "$MODE_ALLSKY" = "both"; then
      if test $CORRELATION != 0.000; then
          SIMSTRING="Simulated                   r = $CORRELATION"
      else
          SIMSTRING="Simulated"
      fi
      if test $AZI_SCALE == 0.10; then
          convert -fill yellow -annotate +19+447 "Simulated"  -pointsize 16 allsky_cyl_$ILOC.png allsky_cyl_$ILOC.png
      elif test $AZI_SCALE == 0.20; then
          convert -fill yellow -annotate +19+179 "Simulated"  -pointsize 14 allsky_cyl_$ILOC.png allsky_cyl_$ILOC.png
      elif test $AZI_SCALE == 0.25; then
          convert -fill yellow -annotate +15+$YDISP  "$SIMSTRING"  -pointsize $POINT allsky_cyl_$ILOC.png allsky_cyl_$ILOC.png
      elif test $AZI_SCALE == 0.50; then
#         convert -fill yellow -annotate +15+$YDISP  "Simulated"  -pointsize 12 allsky_cyl_$ILOC.png allsky_cyl_$ILOC.png
          echo " "
      else
          convert -fill yellow -annotate +15+20  "Simulated"  -pointsize 20 allsky_cyl_$ILOC.png allsky_cyl_$ILOC.png
      fi
  fi
  
# Annotate Time and GHI
  
  if test -n "$ANNOTATE_TIME"; then
      ATIME=$ANNOTATE_TIME
  else
      ATIME=`head -2 label.$ILOC | tail -1`
  fi
  ATIME=`echo $ATIME | sed 's/^[ \t]*//'` # remove leading spaces
  GHIUNITS=W/m^2
  GHI=`head -6 label2.$ILOC | tail -1`$GHIUNITS
  GHI=`echo $GHI | sed 's/^[ \t]*//'`     # remove leading spaces
  echo "Annotate Time $ATIME and GHI $GHI"
  if test "$MODE_ALLSKY" = "polar" || test "$MODE_ALLSKY" = "both"; then
   if test "$IMGGEOM" = "511x511"; then
    echo 'convert -fill white -annotate +360+503  "$ATIME" -pointsize $POINTP allsky_polar_$ILOC.png allsky_polar_$ILOC.png'
          convert -fill white -annotate +360+503  "$ATIME" -pointsize $POINTP allsky_polar_$ILOC.png allsky_polar_$ILOC.png
   elif test "$IMGGEOM" = "1023x1023"; then
    convert -fill white -annotate +770+1003 "$ATIME" -pointsize $POINTP allsky_polar_$ILOC.png allsky_polar_$ILOC.png
   elif test "$IMGGEOM" = "1535x1535"; then
    convert -fill white -annotate +1230+1527 "$ATIME" -pointsize $POINTP allsky_polar_$ILOC.png allsky_polar_$ILOC.png
   elif test "$IMGGEOM" = "2047x2047"; then
    convert -fill white -annotate +1440+2039 "$ATIME" -pointsize $POINTP allsky_polar_$ILOC.png allsky_polar_$ILOC.png
   elif test "$IMGGEOM" = "2559x2559"; then
    convert -fill white -annotate +2300+2531 "$ATIME" -pointsize $POINTP allsky_polar_$ILOC.png allsky_polar_$ILOC.png
   else
    convert -fill white -annotate +1080+1527 "$ATIME" -pointsize $POINTP allsky_polar_$ILOC.png allsky_polar_$ILOC.png
   fi
   echo "POINTP is $POINTP"
  fi

  if test "$MODE_ALLSKY" = "cyl" || test "$MODE_ALLSKY" = "both"; then
    if test $AZI_SCALE == 0.10; then
      convert -fill yellow -annotate +525+447 "$ATIME"  -pointsize 16 allsky_cyl_$ILOC.png allsky_cyl_$ILOC.png
    elif test $AZI_SCALE == 0.20; then
      convert -fill yellow -annotate +918+179 "$ATIME" -pointsize 14 allsky_cyl_$ILOC.png allsky_cyl_$ILOC.png
    elif test $AZI_SCALE == 0.25; then
      XDISP1=480
      XDISP2=815
      if test "$GHI" != ""; then
        echo  convert -fill yellow -annotate  +$XDISP1+$YDISP "$GHI" -annotate +$XDISP2+$YDISP "$ATIME" -pointsize $POINT allsky_cyl_$ILOC.png allsky_cyl_$ILOC.png
              convert -fill yellow -annotate  +$XDISP1+$YDISP "$GHI" -annotate +$XDISP2+$YDISP "$ATIME" -pointsize $POINT allsky_cyl_$ILOC.png allsky_cyl_$ILOC.png
      else
        echo  convert -fill yellow -annotate  +$XDISP2+$YDISP "$ATIME" -pointsize $POINT allsky_cyl_$ILOC.png allsky_cyl_$ILOC.png
              convert -fill yellow -annotate  +$XDISP2+$YDISP "$ATIME" -pointsize $POINT allsky_cyl_$ILOC.png allsky_cyl_$ILOC.png
        echo  convert -fill yellow -annotate  +$XDISP1+$YDISP "$GHI"   -pointsize $POINT allsky_cyl_$ILOC.png allsky_cyl_$ILOC.png
              convert -fill yellow -annotate  +$XDISP1+$YDISP "$GHI"   -pointsize $POINT allsky_cyl_$ILOC.png allsky_cyl_$ILOC.png
      fi
    elif test $AZI_SCALE == 0.50; then
#     convert -fill yellow -annotate +1550+20 "$ATIME" -pointsize 20 allsky_cyl_$ILOC.png allsky_cyl_$ILOC.png
      if test "$GHI" != ""; then
        XDISP1=900
        XDISP2=1550
        convert -fill yellow -annotate  +$XDISP1+$YDISP "$GHI" -annotate +$XDISP2+$YDISP "$ATIME" -pointsize $POINT allsky_cyl_$ILOC.png allsky_cyl_$ILOC.png
      else
        echo 'convert -fill yellow -annotate  +$XDISP2+$YDISP "$ATIME" -pointsize $POINT allsky_cyl_$ILOC.png allsky_cyl_$ILOC.png'
              convert -fill yellow -annotate  +$XDISP2+$YDISP "$ATIME" -pointsize $POINT allsky_cyl_$ILOC.png allsky_cyl_$ILOC.png
        echo 'convert -fill yellow -annotate  +$XDISP1+$YDISP "$GHI"   -pointsize $POINT allsky_cyl_$ILOC.png allsky_cyl_$ILOC.png'
              convert -fill yellow -annotate  +$XDISP1+$YDISP "$GHI"   -pointsize $POINT allsky_cyl_$ILOC.png allsky_cyl_$ILOC.png
      fi
    elif test "$IMGGEOM" = "5761x2881"; then
      XDISP2=3000
      convert -fill white -annotate +$XDISP2+$YDISP "$ATIME" -pointsize 25 allsky_polar_$ILOC.png allsky_polar_$ILOC.png
    else
      convert -fill yellow -annotate +725+20 "$ATIME"  -pointsize 20 allsky_cyl_$ILOC.png allsky_cyl_$ILOC.png
    fi
  fi
  
# Annotate Lat/Lon
# Try and strip off leading blanks?
# LATLON="39.99 -105.26"
  echo "Annotate Lat/Lon $LATLON"
  HT2=`echo $HT  | sed 's/^[0\t]*//'` # remove leading zeros
  HT2=`echo $HT2 | sed 's/\.//'`m     # remove trailing decimal and add 'm'
  LATLON2="$LATLON $HT2"

  if test "$MODE_ALLSKY" = "polar" || test "$MODE_ALLSKY" = "both"; then
   if test "$IMGGEOM" = "511x511"; then
    echo "convert -fill white -annotate +363+20 "$LATLON" -pointsize 18 allsky_polar_$ILOC.png allsky_polar_$ILOC.png"
          convert -fill white -annotate +363+20 "$LATLON" -pointsize 18 allsky_polar_$ILOC.png allsky_polar_$ILOC.png
   elif test "$IMGGEOM" = "1023x1023"; then
    echo "convert -fill white -annotate +875+20 "$LATLON" -pointsize 18 allsky_polar_$ILOC.png allsky_polar_$ILOC.png"
          convert -fill white -annotate +875+20 "$LATLON" -pointsize 18 allsky_polar_$ILOC.png allsky_polar_$ILOC.png
   elif test "$IMGGEOM" = "2559x2559"; then
    echo "convert -fill white -annotate +2300+20 "$LATLON" -pointsize 22 allsky_polar_$ILOC.png allsky_polar_$ILOC.png"
          convert -fill white -annotate +2300+20 "$LATLON" -pointsize 22 allsky_polar_$ILOC.png allsky_polar_$ILOC.png
   else # 1535x1535
    echo "convert -fill white -annotate +1260+20 "$LATLON2" -pointsize $POINTP allsky_polar_$ILOC.png allsky_polar_$ILOC.png"
          convert -fill white -annotate +1260+20 "$LATLON2" -pointsize $POINTP allsky_polar_$ILOC.png allsky_polar_$ILOC.png
   fi
  fi

  if test "$MODE_ALLSKY" = "cyl" || test "$MODE_ALLSKY" = "both"; then
    if test $AZI_SCALE == 0.10; then
      convert -fill yellow -annotate +820+447 "$LATLON" -pointsize 16 allsky_cyl_$ILOC.png allsky_cyl_$ILOC.png
    elif test $AZI_SCALE == 0.20; then
      convert -fill yellow -annotate +1434+179 "$LATLON" -pointsize 14 allsky_cyl_$ILOC.png allsky_cyl_$ILOC.png
    elif test $AZI_SCALE == 0.25; then
      convert -fill yellow -annotate +1202+$YDISP "$LATLON2" -pointsize $POINT allsky_cyl_$ILOC.png allsky_cyl_$ILOC.png
    elif test $AZI_SCALE == 0.50; then
      convert -fill yellow -annotate +1840+20 "$LATLON" -pointsize 20 allsky_cyl_$ILOC.png allsky_cyl_$ILOC.png
    else
      convert -fill yellow -annotate +920+20  "$LATLON" -pointsize 20 allsky_cyl_$ILOC.png allsky_cyl_$ILOC.png
    fi
  fi

# Annotate Field
  if test "$MODE_ALLSKY" = "polar" || test "$MODE_ALLSKY" = "both"; then
    if test "$GHI" != ""; then
      LL=$GHI
      LLPOINT=$POINTP
    else
      LL=All-sky
      LLPOINT=18
    fi

    if test "$IMGGEOM" = "511x511"; then
      convert -fill white -annotate +20+500  "$LL"   -pointsize $POINTP allsky_polar_$ILOC.png allsky_polar_$ILOC.png
    elif test "$IMGGEOM" = "1023x1023"; then
      convert -fill white -annotate +20+1003 "$LL"   -pointsize $POINTP allsky_polar_$ILOC.png allsky_polar_$ILOC.png
    elif test "$IMGGEOM" = "1535x1535"; then
      convert -fill white -annotate +20+1524 "$LL"   -pointsize $POINTP allsky_polar_$ILOC.png allsky_polar_$ILOC.png
    elif test "$IMGGEOM" = "2047x2047"; then
      convert -fill white -annotate +20+1524 "$LL"   -pointsize 20 allsky_polar_$ILOC.png allsky_polar_$ILOC.png
    elif test "$IMGGEOM" = "2559x2559"; then
      convert -fill white -annotate +20+2531 "$LL"   -pointsize 20 allsky_polar_$ILOC.png allsky_polar_$ILOC.png
    else
      convert -fill white -annotate +20+1524 "$LL"   -pointsize 20 allsky_polar_$ILOC.png allsky_polar_$ILOC.png
    fi
  fi
  
# Annotate Directions
  echo "Annotate Directions"
  DIRCYL=South
  DIRCYLL=East
  DIRCYLR=West

  if test "$RESOLUTION_POLAR" = "180pr"; then
      DIR1=SE
      DIR2=NE
      DIR3=NW
      DIR4=SW
  fi
  if test "$RESOLUTION_POLAR" = "360pr"; then # default           
      DIR1=NW
      DIR2=SW
      DIR3=SE
      DIR4=NE
  fi
  if test "$RESOLUTION_POLAR" = "360p"; then # flip left/right
      DIR1=NE
      DIR2=SE
      DIR3=SW
      DIR4=NW
  fi
  if test "$RESOLUTION_POLAR" = "180p"; then # rotate and flip left/right
      DIR1=SW
      DIR2=NW
      DIR3=NE
      DIR4=SE
  fi
  if test "$RESOLUTION_CYL" = "360c"; then # roll horizontally by half the image
      DIRCYL=North
      DIRCYLL=West
      DIRCYLR=East
  fi

# if test "$MODE_ALLSKY" = "polar" || test "$MODE_ALLSKY" = "both"; then
#  if test "$IMGGEOM" = "511x511"; then
#   convert -fill white -annotate +55+60     "$DIR1"            -pointsize 20 allsky_polar_$ILOC.png allsky_polar_$ILOC.png
#   convert -fill white -annotate +40+450    "$DIR2"            -pointsize 20 allsky_polar_$ILOC.png allsky_polar_$ILOC.png
#   convert -fill white -annotate +440+450   "$DIR3"            -pointsize 20 allsky_polar_$ILOC.png allsky_polar_$ILOC.png
#   convert -fill white -annotate +435+60    "$DIR4"            -pointsize 20 allsky_polar_$ILOC.png allsky_polar_$ILOC.png
#  fi
#  fi

  POINTCYL=16

  if test "$MODE_ALLSKY" = "cyl" || test "$MODE_ALLSKY" = "both"; then
    if test $AZI_SCALE == 0.10; then
      convert -fill orange -annotate  +1050+447   "East"          -pointsize 16 allsky_cyl_$ILOC.png allsky_cyl_$ILOC.png
    elif test $AZI_SCALE == 0.20; then
      convert -fill orange -annotate  +776+179  "$DIRCYL"          -pointsize 14 allsky_cyl_$ILOC.png allsky_cyl_$ILOC.png
    elif test $AZI_SCALE == 0.25; then
      convert -fill orange -annotate  +693+$YDISP   "$DIRCYL"          -pointsize $POINTCYL allsky_cyl_$ILOC.png allsky_cyl_$ILOC.png
      convert -fill orange -annotate  +339+$YDISP   "$DIRCYLL"         -pointsize $POINTCYL allsky_cyl_$ILOC.png allsky_cyl_$ILOC.png
      convert -fill orange -annotate  +1055+$YDISP  "$DIRCYLR"         -pointsize $POINTCYL allsky_cyl_$ILOC.png allsky_cyl_$ILOC.png
    elif test $AZI_SCALE == 0.50; then
      convert -fill orange -annotate +1040+20   "$DIRCYL"          -pointsize 20 allsky_cyl_$ILOC.png allsky_cyl_$ILOC.png
    else
      convert -fill orange -annotate +520+20    "$DIRCYL"          -pointsize 20 allsky_cyl_$ILOC.png allsky_cyl_$ILOC.png
    fi
  fi
  
