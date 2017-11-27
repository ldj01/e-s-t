; Kelly G. Laraby
; 12 March 2015
;
;
; NAME:  lst_cfsr_step3
;
; PURPOSE:  IDL PROCEDURE
;   Generate transmission, upwelled radiance, and downwelled radiance at each pixel
;
; CALL SEQUENCE: lst_merra_step3,  home, $            ;directory containing programs and supporting files
;                                  directory, $       ;directory containing landsat metadata and location for results
;                                  imageBase, $       ;basename of landsat files
;                                  numPoints, $       ;number of CFSR points within Landsat scene
;                                  numHeights, $      ;number of height at each CFSR point
;                                  landsatSamples, $  ;number of columns in thermal landsat scene
;                                  landsatLines, $    ;number of rows in thermal landsat scene
;                                  pixelSize, $       ;size of pixel in thermal landsat scene (in meters)
;                                  zone, $            ;utm zone of landsat scene (for utm corrdinates)
;                                  UL_EAST, $         ;easting coordinate of upper left corner
;                                  UL_NORTH, $        ;northing coordiante of upper left corner
;                                  UR_EAST, $         ;easting coordinate of upper right corner
;                                  UR_NORTH, $        ;northing coordinate of upper right corner
;                                  LL_EAST, $         ;easting coordinate of lower left corner
;                                  LL_NORTH, $        ;northing coordinate of lower left corner
;                                  LR_EAST, $         ;easting coordinate of lower right corner
;                                  LR_NORTH,$         ;northing coordinate of lower right corner
;                                  UL_LAT, $          ;latitude of upper left corner
;                                  UL_LON, $          ;longitude of upper left corner
;                                  UR_LAT, $          ;latitude of upper right corner
;                                  UR_LON, $          ;longitude of upper right corner
;                                  LL_LAT, $          ;latitude of lower left corner
;                                  LL_LON, $          ;longitude of lower left corner
;                                  LR_LAT, $          ;latitude of lower right corner
;                                  LR_LON, $          ;longitude of lower right corner
;                                  LMAX6, $           ;calibration coefficient for Landsat Thermal band
;                                  LMIN6, $           ;calibration coefficient for Landsat Thermal band
;                                  QCALMAX6, $        ;calibration coefficient for Landsat Thermal band
;                                  QCALMIN6, $        ;calibration coefficient for Landsat Thermal band
;                                  demFile, $         ;file path containing DEM for current scene
;                                  whichLandsat       ;which Landsat sensors was this image captured with
;
;  RESTRICTIONS:
;    Already called lst_cfsr_step1.pro, commandList, parsed tape6 files, and lst_cfsr_step2.pro
;
;  REQUIRED PROGRAMS AND FILES (in home directory):   atmosphericParameters.txt
;                                                     CONVERT_LL_UTM.pro
;                                                     DISTANCE_IN_UTM.pro
;                                                     INTERP_TO_HEIGHT.pro
;                                                     INTERP_TO_LOCATION.pro
;                                                     coordinates.txt
;
;

PRO lst_cfsr_step3,  home, $
  directory, $
  imageBase, $
  numPoints, $
  numHeights, $
  landsatSamples, $
  landsatLines, $
  pixelSize, $
  zone, $
  UL_EAST, $
  UL_NORTH, $
  UR_EAST, $
  UR_NORTH, $
  LL_EAST, $
  LL_NORTH, $
  LR_EAST, $
  LR_NORTH, $
  UL_LAT, $
  UL_LON, $
  UR_LAT, $
  UR_LON, $
  LL_LAT, $
  LL_LON, $
  LR_LAT, $
  LR_LON, $
  LMAX6, $
  LMIN6, $
  QCALMAX6, $
  QCALMIN6, $
  demFile, $
  whichLandsat

  ;convert inputs
  home = STRING(home)
  directory = STRING(directory)
  imageBase = STRING(imageBase)
  numPoints = FIX(numPoints)
  numHeights = FIX(numHeights)
  landsatSamples = FLOAT(landsatSamples)
  landsatLines = FLOAT(landsatLines)
  pixelSize = FLOAT(pixelSize)
  zone = FLOAT(zone)
  UL_EAST = FLOAT(UL_EAST)
  UL_NORTH = FLOAT(UL_NORTH)
  UR_EAST = FLOAT(UR_EAST)
  UR_NORTH = FLOAT(UR_NORTH)
  LL_EAST = FLOAT(LL_EAST)
  LL_NORTH = FLOAT(LL_NORTH)
  LR_EAST = FLOAT(LR_EAST)
  LR_NORTH = FLOAT(LR_NORTH)
  UL_LAT = FLOAT(UL_LAT)
  UL_LON = FLOAT(UL_LON)
  UR_LAT = FLOAT(UR_LAT)
  UR_LON = FLOAT(UR_LON)
  LL_LAT = FLOAT(LL_LAT)
  LL_LON = FLOAT(LL_LON)
  LR_LAT = FLOAT(LR_LAT)
  LR_LON = FLOAT(LR_LON)
  LMAX6 = FLOAT(LMAX6)
  LMIN6 = FLOAT(LMIN6)
  QCALMAX6 = FLOAT(QCALMAX6)
  QCALMIN6 = FLOAT(QCALMIN6)
  demFile = STRING(demFile)
  whichLandsat = FIX(whichLandsat)

  ;define file containg atmospheric parameters
  atmFile = home+directory+'atmosphericParameters.txt'

  ;define file containing landsat image band 6
  CASE whichLandsat OF
    5: landsatThermalFile = home+directory+imageBase+'_B6.TIF'
    7: landsatThermalFile = home+directory+imageBase+'_B6_VCID_1.TIF'
    10: landsatThermalFile = home+directory+imageBase+'_B10.TIF'
    11: landsatThermalFile = home+directory+imageBase+'_B11.TIF'
  ENDCASE

  ;read in landsat image
  landsatBand6 = READ_TIFF(landsatThermalFile, GEOTIFF=landsatGeotiff)
  zeros = WHERE(landsatBand6 EQ 0)
  mask = MAKE_ARRAY(landsatSamples, landsatLines, VALUE = 1)
  mask[zeros] = 0
  ;convert digital counts to radiance for landsat band 6 [ W m^(-2) sr^(-1) mu^(-1) ]
  landsatThermal = [((LMAX6 - LMIN6)/(QCALMAX6 - QCALMIN6))*(landsatBand6 - QCALMIN6)+LMIN6]
  IF whichLandsat EQ 5 THEN landsatThermal = landsatThermal+0.044

  landsatThermal = landsatThermal*mask

  ;read in DEM for current image (values in meters)
  PRINT, "looking for dem in:", demFile
  dem = READ_TIFF(demFile, GEOTIFF=demGeotiff)

  ;easting and northing values of the upper left corner of the image including the fill space.
  UL_demEasting = demGeotiff.MODELTIEPOINTTAG[3]
  UL_demNorthing = demGeotiff.MODELTIEPOINTTAG[4]

  ;determine northern or southern hemisphere. '6' = northern hemisphere, '7' = southern hermisphere.
  demHemi = STRMID(STRCOMPRESS(STRING(demGEOTIFF.PROJECTEDCSTYPEGEOKEY)),3,1)
  demHemi = FIX(demHemi)

  IF demHemi EQ 7 THEN UL_demNorthing = UL_demNorthing - 10000000D

  ;determine dem UTM zone
  ;   demUTM = STRMID(STRCOMPRESS(STRING(demGeotiff.PROJECTEDCSTYPEGEOKEY)),4,2)
  ;   demZone = demUTM

  ;determine pixel size of dem
  demPixelX = demGeotiff.MODELPIXELSCALETAG[0]
  demPixelY = demGeotiff.MODELPIXELSCALETAG[1]

  ;determine offset in integer number of pixels between dem and landsat scene
  offsetEasting = ROUND((UL_EAST-UL_demEasting)/30D)
  offsetNorthing = ROUND(ABS(UL_demNorthing-UL_NORTH)/30D)

  ;read in atmospheric parameters
  OPENR, 20, atmFile
  atm = MAKE_ARRAY(6, numPoints*numHeights, /FLOAT)
  READF, 20, atm
  CLOSE, 20
  FREE_LUN, 20

  ;determine if landsat is in the northern or southern hemisphere. '6' = northern hemisphere, '7' = southern hermisphere.
  IF UL_LAT GE 0 THEN landsatHemi = 6 ELSE landsatHemi = 7

  ;define sampling for latitude and longitude of merra grid
  deltalambda = 0.5
  deltaphi = 0.5

  ;
  ;define lon grid for CFSR data
  ;
  i = INDGEN(720)+1
  lonarray = -180 + deltalambda*(i-1.0)
  lon = MAKE_ARRAY(720,361, /DOUBLE)
  eye = MAKE_ARRAY(720,361, /DOUBLE)
  FOR g = 0, 360 DO BEGIN
    lon[*,g] = lonarray
    eye[*,g] = i
  ENDFOR

  ;
  ;define lat grid for CFSR data
  ;
  j = INDGEN(361)+1
  latarray = -90 + deltaphi*(j-1.0)
  lat = MAKE_ARRAY(720,361, /DOUBLE)
  jay = MAKE_ARRAY(720,361, /DOUBLE)
  FOR g = 0, 719 DO BEGIN
    lat[g,*] = latarray
    jay[g,*] = j
  ENDFOR

  ;
  ; expand range to include NARR points outside image for edge pixels
  ;
  UL_LAT = UL_LAT + 0.6 ;1.3
  UL_LON = UL_LON - 0.6 ;1.3
  LR_LAT = LR_LAT - 0.6 ;1.3
  LR_LON = LR_LON + 0.6 ;1.3

  ;determine what points in the MERRA dataset fall within the Landsat image using logical operators
  ;lessThanLat and greaterThanLat are values where the MERRA values are less than or greater than the edges of the Landsat
  ;corners values respectively
  ;pixels that are true in both fall within the Landsat scene
  ;the same thing is done with longitude values
  lessThanLat = DOUBLE(lat LT UL_LAT)
  greaterThanLat = DOUBLE(lat GT LR_LAT)
  keepLat = lessThanLat*greaterThanLat
  greaterThanLon = DOUBLE(lon GT UL_LON)
  lessThanLon = DOUBLE(lon LT LR_LON)
  keepLon = greaterThanLon*lessThanLon

  ;values that are true in both keepLat and keepLon fall within the Landsat image
  ;convert indices into (x,y) values in the MERRA dataset
  ;Because the Landsat is in easting/northing originally and the MERRA data is lat/lon, this is not a perfect square
  keep = keepLat * keepLon
  inLandsat = WHERE(keep NE 0)

  ;determine indices to pull out rectangle of MERRA points
  iindices = [MIN(eye[inLandsat])-1,MAX(eye[inLandsat])-1]
  jindices = [MIN(jay[inLandsat])-1,MAX(jay[inLandsat])-1]

  PRINT, 'size of lat = ',SIZE(lat)
  PRINT, 'i ind = ', iindices[0], iindices[1]
  PRINT, 'j ind = ', jindices[0], jindices[1]
  ;extract coordinates within this rectangle
  igrid = eye[iindices[0]:iindices[1],jindices[0]:jindices[1]]
  jgrid = jay[iindices[0]:iindices[1],jindices[0]:jindices[1]]
  latgrid = lat[iindices[0]:iindices[1],jindices[0]:jindices[1]]
  longrid = lon[iindices[0]:iindices[1],jindices[0]:jindices[1]]

  ivector = REFORM(igrid,1,N_ELEMENTS(igrid))
  jvector = REFORM(jgrid,1,N_ELEMENTS(jgrid))

  inLandsatArray = [ivector-1,jvector-1]


  size = SIZE(igrid, /DIMENSIONS)
  cfsrSamples = size[0]
  cfsrLines = size[1]

  cfsrUTM = CONVERT_LL_UTM(latgrid, longrid, zone)

  eastvector = cfsrUTM[1,*]
  northvector = cfsrUTM[2,*]
  eastGrid = REFORM(cfsrUTM[1,*],cfsrSamples,cfsrLines)
  northGrid = REFORM(cfsrUTM[2,*],cfsrSamples,cfsrLines)

  ;initialize arrays to save results into
  transmission = MAKE_ARRAY(landsatSamples, landsatLines, /FLOAT)
  upwelled = MAKE_ARRAY(landsatSamples, landsatLines, /FLOAT)
  downwelled = MAKE_ARRAY(landsatSamples, landsatLines, /FLOAT)
  elevation = MAKE_ARRAY(landsatSamples, landsatLines, /FLOAT)
  upperLeftArray = MAKE_ARRAY(2, landsatSamples, landsatLines, /FLOAT)

  ;stoph = !VALUES.D_INFINITY
  ;stopk = !VALUES.D_INFINITY

  ;iterature through all rows in landsat scene
  FOR k = 0, landsatLines-1 DO BEGIN

    ;set first in line to true for each row
    firstInLine = 1

    ;iterature through all columns in landsat scene
    FOR h = 0, landsatSamples - 1 DO BEGIN

      IF landsatThermal[h,k] NE 0 THEN BEGIN

        ;determine UTM coordinates of current pixel
        currentEasting = UL_EAST + h*pixelSize
        currentNorthing = UL_NORTH - k*pixelSize

        IF firstInLine THEN BEGIN

          ;compute distance between current pixel and each CFSR points in UTM coordinates
          distances = MAKE_ARRAY(numPoints, /FLOAT)
          FOR g = 0, numPoints - 1 DO BEGIN
            distances[g] = DISTANCE_IN_UTM(eastvector[g], northvector[g], currentEasting, currentNorthing)
          ENDFOR

          ;find indices of 6 closest points
          closest = SORT(distances)
          closest = closest[0:6]

          ;extract closest points and determine if they are above or below
          eastingsNear = eastvector[closest]
          northingsNear = northvector[closest]
          above = WHERE(northingsNear GT currentNorthing)
          below = WHERE(northingsNear LE currentNorthing)

          ;find upper left point of quad
          upperLeft = [MIN(ivector[closest[below]]), MIN(jvector[closest[below]])]
          ;find indices in grid of NARR points
          upperLeftIndices = ARRAY_INDICES(igrid,WHERE(igrid EQ upperLeft[0] AND jgrid EQ upperLeft[1]))

          ;extract UTM coordinates of four points to be interpolated and build array
          UTMa = [eastGrid[upperLeftIndices[0],upperLeftIndices[1]],northGrid[upperLeftIndices[0],upperLeftIndices[1]]]
          UTMb = [eastGrid[upperleftIndices[0]+1,upperLeftIndices[1]],northGrid[upperleftIndices[0]+1,upperLeftIndices[1]]]
          UTMc = [eastGrid[upperleftIndices[0],upperLeftIndices[1]+1],northGrid[upperleftIndices[0],upperLeftIndices[1]+1]]
          UTMd = [eastGrid[upperleftIndices[0]+1,upperLeftIndices[1]+1],northGrid[upperleftIndices[0]+1,upperLeftIndices[1]+1]]
          coordinates = [[UTMa],[UTMb],[UTMc],[UTMd]]

          ;determine index of four points in order to pull from atmospheric parameter file
          indexA = WHERE(eastGrid EQ UTMa[0] AND northGrid EQ UTMa[1])
          indexB = WHERE(eastGrid EQ UTMb[0] AND northGrid EQ UTMb[1])
          indexC = WHERE(eastGrid EQ UTMc[0] AND northGrid EQ UTMc[1])
          indexD = WHERE(eastGrid EQ UTMd[0] AND northGrid EQ UTMd[1])
          indices = [indexA, indexB, indexC, indexD]

          ;set firstInLine variable to false
          firstInLine = 0

        ENDIF ELSE BEGIN

          ;given indices of previous pixel, there are six possible quads to move into
          ;check 6 distances to determine new upperleft corner

          stayRight = DISTANCE_IN_UTM(eastGrid[upperLeftIndices[0],upperLeftIndices[1]], $
            northGrid[upperLeftIndices[0],upperLeftIndices[1]], $
            currentEasting, currentNorthing)

          IF upperLeftIndices[0]+2 LT cfsrSamples THEN BEGIN
            moveRight = DISTANCE_IN_UTM(eastGrid[upperLeftIndices[0]+2,upperLeftIndices[1]], $
              northGrid[upperLeftIndices[0]+2,upperLeftIndices[1]], $
              currentEasting, currentNorthing)
          ENDIF ELSE BEGIN
            moveRight = !VALUES.D_INFINITY
          ENDELSE

          stayUp = DISTANCE_IN_UTM(eastGrid[upperLeftIndices[0]+1,upperLeftIndices[1]+1], $
            northGrid[upperLeftIndices[0]+1,upperLeftIndices[1]+1], $
            currentEasting, currentNorthing)
          moveUp = DISTANCE_IN_UTM(eastGrid[upperLeftIndices[0]+1,upperLeftIndices[1]-1], $
            northGrid[upperLeftIndices[0]+1,upperLeftIndices[1]-1], $
            currentEasting, currentNorthing)

          stayDown = DISTANCE_IN_UTM(eastGrid[upperLeftIndices[0]+1,upperLeftIndices[1]], $
            northGrid[upperLeftIndices[0]+1,upperLeftIndices[1]], $
            currentEasting, currentNorthing)
          IF upperLeftIndices[1]+2 LT cfsrLines THEN BEGIN
            moveDown = DISTANCE_IN_UTM(eastGrid[upperLeftIndices[0]+1,upperLeftIndices[1]+2], $
              northGrid[upperLeftIndices[0]+1,upperLeftIndices[1]+2], $
              currentEasting, currentNorthing)
          ENDIF ELSE BEGIN
            moveDown = !VALUES.D_INFINITY
          ENDELSE

          IF moveRight LT stayRight THEN upperLeftIndices = [upperLeftIndices[0]+1,upperLeftIndices[1]]
          IF moveUp LT stayUp THEN upperLeftIndices = [upperLeftIndices[0],upperLeftIndices[1]-1]
          IF moveDown LT stayDown THEN upperLeftIndices = [upperLeftIndices[0],upperLeftIndices[1]+1]
          upperLeft = [igrid[upperLeftIndices[0],upperLeftIndices[1]],jgrid[upperLeftIndices[0],upperLeftIndices[1]]]

          ;extract UTM coordinates of four points to be interpolated and build array
          UTMa = [eastGrid[upperLeftIndices[0],upperLeftIndices[1]],northGrid[upperLeftIndices[0],upperLeftIndices[1]]]
          UTMb = [eastGrid[upperleftIndices[0]+1,upperLeftIndices[1]],northGrid[upperleftIndices[0]+1,upperLeftIndices[1]]]
          UTMc = [eastGrid[upperleftIndices[0],upperLeftIndices[1]+1],northGrid[upperleftIndices[0],upperLeftIndices[1]+1]]
          UTMd = [eastGrid[upperleftIndices[0]+1,upperLeftIndices[1]+1],northGrid[upperleftIndices[0]+1,upperLeftIndices[1]+1]]
          coordinates = [[UTMa],[UTMb],[UTMc],[UTMd]]

          ;determine index of four points in order to pull from atmospheric parameter file
          indexA = WHERE(eastGrid EQ UTMa[0] AND northGrid EQ UTMa[1])
          indexB = WHERE(eastGrid EQ UTMb[0] AND northGrid EQ UTMb[1])
          indexC = WHERE(eastGrid EQ UTMc[0] AND northGrid EQ UTMc[1])
          indexD = WHERE(eastGrid EQ UTMd[0] AND northGrid EQ UTMd[1])
          indices = [indexA, indexB, indexC, indexD]

        ENDELSE

        ;print, h, k
        ;IF h EQ 1000 AND k EQ 1000 THEN stop
        ;IF h EQ 2000 AND k EQ 2000 THEN stop
        ;IF h EQ 3000 AND k EQ 3000 THEN stop
        ;IF h EQ 4000 AND k EQ 4000 THEN stop
        ;IF h EQ 5000 AND k EQ 5000 THEN stop
        ;IF h EQ 6000 AND k EQ 6000 THEN stop

        upperLeftArray[0,h,k] = upperLeftIndices[0]
        upperLeftArray[1,h,k] = upperLeftIndices[1]

        ;get elevation of current pixel from dem
        height = dem[h+offsetEasting, k+offsetNorthing]
        ;convert height from m to km
        height = height/1000D

        ;create matrix to hold atmospheric parameters interpolated to appropriate height
        atHeight = MAKE_ARRAY(3, N_ELEMENTS(indices), /FLOAT)

        ;interpolate three parameters to that height at each of the four closest points
        FOR g = 0, N_ELEMENTS(indices)-1 DO BEGIN

          ;extract atmospheric parameters for all heights at current location
          currentIndex = indices[g]*numHeights
          currentLocation = atm[*,currentIndex:currentIndex+numHeights-1]

          ;separate height and atmospheric data
          currentHeight = currentLocation[2,*]
          currentAtm = currentLocation[3:5,*]

          ;interpolate three atmospheric parameters to current height
          atHeight[0:2,g] = INTERP_TO_HEIGHT(currentHeight, currentAtm, height)

        ENDFOR

        ;interpolate parameters at appropriate height to location of current pixel
        parameters = INTERP_TO_LOCATION(coordinates, atHeight, currentEasting, currentNorthing)

        ;save atmospheric parameters to full arrays
        transmission[h,k] = parameters[0]
        upwelled[h,k] = parameters[1]
        downwelled[h,k] = parameters[2]
        elevation[h,k] = height

      ENDIF

      ;end for loop iterating through samples
    ENDFOR

    print, k

    ;end for loop iterating through lines
  ENDFOR

  ;convert radiances to W m^(-2) sr^(-1)
  upwelled = upwelled*(100)^2
  downwelled = downwelled*(100)^2
  elevation = elevation*(1000)

  ;save results in save variable
  ;SAVE, landsatThermal, elevation, transmission, upwelled, downwelled, FILENAME = home+directory+'variables.sav'

  ;write results to float precision tiff
  ;   results = MAKE_ARRAY(5, landsatSamples, landsatLines, /DOUBLE)
  ;   results[0,*,*] = landsatThermal
  ;   results[1,*,*] = elevation
  ;   results[2,*,*] = transmission
  ;   results[3,*,*] = upwelled
  ;   results[4,*,*] = downwelled

  SAVE, upperLeftArray, FILENAME=home+directory+'upperLeft.sav'

  WRITE_TIFF, home+directory+'upperLeft.tif', upperLeftArray, GEOTIFF=landsatGeotiff, /FLOAT

  WRITE_TIFF, home+directory+imagebase+'_LSTresults.tif', $
    TRANSPOSE([[[landsatThermal]],[[elevation]],[[transmission]],[[upwelled]],[[downwelled]]],[2,0,1]), $
    GEOTIFF=landsatGeotiff, $
    /FLOAT

END
