/*****************************************************************************
!File: input.c
*****************************************************************************/

#include "input.h"

/* Functions */
Input_t *OpenInput(Espa_internal_meta_t *metadata)
/* 
!C******************************************************************************

!Description: 'OpenInput' sets up the 'input' data structure, opens the
 input raw binary files for read access.
 
!Input Parameters:
 metadata     'Espa_internal_meta_t' data structure with XML info

!Output Parameters:
 (returns)      'input' data structure or NULL when an error occurs

!Team Unique Header:

!END****************************************************************************
*/
{
  Input_t *this = NULL;
  char *error_string = NULL;

  /* Create the Input data structure */
  this = (Input_t *)malloc(sizeof(Input_t));
  if (this == NULL) 
    RETURN_ERROR("allocating Input data structure", "OpenInput", NULL);

  /* Initialize and get input from header file */
  if (!GetXMLInput (this, metadata)) {
    free(this);
    this = NULL;
    RETURN_ERROR("getting input from header file", "OpenInput", NULL);
  }

  /* Open files for access */
  if (this->file_type == INPUT_TYPE_BINARY) {
    if ( this->nband_th == 1 || this->nband_th == 2 ) {
      this->fp_bin_th = fopen(this->file_name_th, "r");
      if (this->fp_bin_th == NULL) 
        error_string = "opening thermal binary file";
      else
        this->open_th = true;
    }
  } else 
    error_string = "invalid file type";

  if (error_string != NULL) {
    free(this->file_name_th);
    this->file_name_th = NULL;
    if ( this->file_type == INPUT_TYPE_BINARY )
      fclose(this->fp_bin_th);  
    this->open_th = false;
    free(this);
    this = NULL;
    RETURN_ERROR(error_string, "OpenInput", NULL);
  }

  return this;
}

bool GetInputLineTh(Input_t *this, int iline, unsigned char *line) 
{
  long loc;
  void *buf_void = NULL;

  if (this == NULL) 
    RETURN_ERROR("invalid input structure", "GetInputLine", false);
  if ( this->nband_th < 1 ) 
    RETURN_ERROR("no thermal input band", "GetInputLine", false);
  if (iline < 0  ||  iline >= this->size_th.l) 
    RETURN_ERROR("line index out of range", "GetInputLine", false);
  if (!this->open_th)
    RETURN_ERROR("band not open", "GetInputLine", false);

  buf_void = (void *)line;
  if (this->file_type == INPUT_TYPE_BINARY) {
    loc = (long) (iline * this->size_th.s * sizeof(uint8));
    if (fseek(this->fp_bin_th, loc, SEEK_SET))
      RETURN_ERROR("error seeking line (binary)", "GetInputLine", false);
    if (fread(buf_void, sizeof(uint8), (size_t)this->size_th.s, 
              this->fp_bin_th) != (size_t)this->size_th.s)
      RETURN_ERROR("error reading line (binary)", "GetInputLine", false);
  }

  return true;
}


bool CloseInput(Input_t *this)
/* 
!C******************************************************************************

!Description: 'CloseInput' ends SDS access and closes the input file.
 
!Input Parameters:
 this           'input' data structure

!Output Parameters:
 this           'input' data structure; the following fields are modified:
                   open
 (returns)      status:
                  'true' = okay
                  'false' = error return

!Team Unique Header:

!END****************************************************************************
*/
{
  if (this == NULL) 
    RETURN_ERROR("invalid input structure", "CloseInput", false);

  /*** now close the thermal file ***/
  if (this->open_th) 
  {
    if (this->file_type == INPUT_TYPE_BINARY)
      fclose(this->fp_bin_th);
    this->open_th = false;
  }

  return true;
}


bool FreeInput(Input_t *this)
/* 
!C******************************************************************************

!Description: 'FreeInput' frees the 'input' data structure memory.
 
!Input Parameters:
 this           'input' data structure

!Output Parameters:
 (returns)      status:
                  'true' = okay (always returned)

!Team Unique Header:

!END****************************************************************************
*/
{
    free(this->file_name_th);
    this->file_name_th = NULL;

    free(this);
    this = NULL;

  return true;
}

#define DATE_STRING_LEN (50)
#define TIME_STRING_LEN (50)

bool GetXMLInput(Input_t *this, Espa_internal_meta_t *metadata)
/* 
!C******************************************************************************

!Description: 'GetXMLInput' pulls input values from the XML structure.
 
!Input Parameters:
 this         'Input_t' data structure to be populated
 metadata     'Espa_internal_meta_t' data structure with XML info

!Output Parameters:
 (returns)      status:
                  'true' = okay (always returned)
                  'false' = error getting metadata from the XML file

!Team Unique Header:

! Design Notes:
  1. This replaces the previous GetInputMeta so the input values are pulled
     from the XML file instead of the HDF and MTL files.
!END****************************************************************************
*/
{
    char *error_string = NULL;
    char acq_date[DATE_STRING_LEN + 1];
    char prod_date[DATE_STRING_LEN + 1];
    char acq_time[TIME_STRING_LEN + 1];
    char temp[MAX_STR_LEN + 1];
    int th_indx;       /* band index in XML file for the thermal band */
    Espa_global_meta_t *gmeta = &metadata->global; /* pointer to global meta */

    /* Initialize the input fields.  Set file type to binary, since that is
       the ESPA internal format for the input L1G/T products. */
    this->file_type = INPUT_TYPE_BINARY;
    this->meta.sat = SAT_NULL;
    this->meta.inst = INST_NULL;
    this->meta.acq_date.fill = true;
    this->meta.time_fill = true;
    this->meta.prod_date.fill = true;
    this->meta.sun_zen = ANGLE_FILL;
    this->meta.sun_az = ANGLE_FILL;
    this->meta.wrs_sys = (Wrs_t)WRS_FILL;
    this->meta.ipath = -1;
    this->meta.irow = -1;
    this->meta.fill = INPUT_FILL;
    this->size_th.s = this->size_th.l = -1;
    this->nband_th = 0;
    this->open_th = false;
    this->meta.gain_th = GAIN_BIAS_FILL;
    this->meta.bias_th = GAIN_BIAS_FILL;
    this->file_name_th = NULL;
    this->fp_bin_th = NULL;

    /* Pull the appropriate data from the XML file */
    acq_date[0] = acq_time[0] = '\0';
    prod_date[0] = '\0';
    if (!strcmp (gmeta->satellite, "LANDSAT_1"))
        this->meta.sat = SAT_LANDSAT_1;
    else if (!strcmp (gmeta->satellite, "LANDSAT_2"))
        this->meta.sat = SAT_LANDSAT_2;
    else if (!strcmp (gmeta->satellite, "LANDSAT_3"))
        this->meta.sat = SAT_LANDSAT_3;
    else if (!strcmp (gmeta->satellite, "LANDSAT_4"))
        this->meta.sat = SAT_LANDSAT_4;
    else if (!strcmp (gmeta->satellite, "LANDSAT_5"))
        this->meta.sat = SAT_LANDSAT_5;
    else if (!strcmp (gmeta->satellite, "LANDSAT_7"))
        this->meta.sat = SAT_LANDSAT_7;
    else if (!strcmp (gmeta->satellite, "LANDSAT_8"))
        this->meta.sat = SAT_LANDSAT_8;
    else
    {
        sprintf (temp, "invalid satellite; value = %s", gmeta->satellite);
        RETURN_ERROR (temp, "GetXMLInput", true);
    }

    if (!strcmp (gmeta->instrument, "TM"))
        this->meta.inst = INST_TM;
    else if (!strncmp (gmeta->instrument, "ETM", 3))
        this->meta.inst = INST_ETM;
    else if (!strncmp (gmeta->instrument, "OLI_TIRS", 8))
        this->meta.inst = INST_ETM;
    else
    {
        sprintf (temp, "invalid instrument; value = %s", gmeta->instrument);
        RETURN_ERROR (temp, "GetXMLInput", true);
    }

    strcpy (acq_date, gmeta->acquisition_date);
    strcpy (acq_time, gmeta->scene_center_time);
    this->meta.time_fill = false;

    /* Make sure the acquisition time is not too long (i.e. contains too
       many decimal points for the date/time routines).  The time should be
       hh:mm:ss.ssssssZ (see DATE_FORMAT_DATEA_TIME in date.h) which is 16
       characters long.  If the time is longer than that, just chop it off. */
    if (strlen (acq_time) > 16)
        sprintf (&acq_time[15], "Z");

    this->meta.sun_zen = gmeta->solar_zenith;
    if (this->meta.sun_zen < -90.0 || this->meta.sun_zen > 90.0)
    {
        error_string = "solar zenith angle out of range";
        RETURN_ERROR (error_string, "GetXMLInput", true);
    }
    this->meta.sun_zen *= RAD;   /* convert to radians */

    this->meta.sun_az = gmeta->solar_azimuth;
    if (this->meta.sun_az < -360.0 || this->meta.sun_az > 360.0)
    {
        error_string = "solar azimuth angle out of range";
        RETURN_ERROR (error_string, "GetXMLInput", true);
    }
    this->meta.sun_az *= RAD;    /* convert to radians */

    switch (gmeta->wrs_system)
    {
        case 1: this->meta.wrs_sys = WRS_1; break;
        case 2: this->meta.wrs_sys = WRS_2; break;
        default:
            sprintf (temp, "invalid WRS system; value = %d",
                gmeta->wrs_system);
            RETURN_ERROR (temp, "GetXMLInput", true);
    }
    this->meta.ipath = gmeta->wrs_path;
    this->meta.irow = gmeta->wrs_row;

    if (this->meta.inst == INST_TM || this->meta.inst == INST_ETM)
    {
        this->nband_th = 1;  /* number of thermal bands; only use 6L for ETM */
        this->meta.iband_th[0] = 6;
        th_indx = 5;
        this->meta.gain_th = metadata->band[th_indx].toa_gain;
        this->meta.bias_th = metadata->band[th_indx].toa_bias;
    }
    else /* this->meta.inst == INST_OLI_TIRS */
    {
        this->nband_th = 2;     /* number of thermal bands */
        this->meta.iband_th[0] = 10;
        this->meta.iband_th[1] = 11;
        th_indx = 9;
        this->meta.gain_th = metadata->band[th_indx].toa_gain;
        this->meta.bias_th = metadata->band[th_indx].toa_bias;
    }

    /* Pull the reflectance info from thermal in the XML file */
    this->size_th.s = metadata->band[th_indx].nsamps;
    this->size_th.l = metadata->band[th_indx].nlines;
    this->meta.pixel_size[0] = metadata->band[th_indx].pixel_size[0];
    this->meta.pixel_size[1] = metadata->band[th_indx].pixel_size[1];
    this->meta.scale_factor_th = metadata->band[th_indx].scale_factor;

    /* Check WRS path/rows */
    if (this->meta.wrs_sys == WRS_1)
    {
        if (this->meta.ipath > 251)
            error_string = "WRS path number out of range";
        else if (this->meta.irow > 248)
            error_string = "WRS row number out of range";
    }
    else if (this->meta.wrs_sys == WRS_2)
    {
        if (this->meta.ipath > 233)
            error_string = "WRS path number out of range";
        else if (this->meta.irow > 248)
            error_string = "WRS row number out of range";
    }
    else
        error_string = "invalid WRS system";

    if (error_string != NULL)
    {
        RETURN_ERROR (error_string, "GetHeaderInput", true);
    }

    /* Check satellite/instrument combination */
    if (this->meta.inst == INST_MSS)
    {
        if (this->meta.sat != SAT_LANDSAT_1 &&
            this->meta.sat != SAT_LANDSAT_2 &&
            this->meta.sat != SAT_LANDSAT_3 &&
            this->meta.sat != SAT_LANDSAT_4 &&
            this->meta.sat != SAT_LANDSAT_5)
            error_string = "invalid insturment/satellite combination";
    }
    else if (this->meta.inst == INST_TM)
    {
        if (this->meta.sat != SAT_LANDSAT_4 &&
            this->meta.sat != SAT_LANDSAT_5)
            error_string = "invalid insturment/satellite combination";
    }
    else if (this->meta.inst == INST_ETM)
    {
        if (this->meta.sat != SAT_LANDSAT_7)
            error_string = "invalid insturment/satellite combination";
    }
    else if (this->meta.inst == INST_OLI_TIRS)
    {
        if (this->meta.sat != SAT_LANDSAT_8)
            error_string = "invalid insturment/satellite combination";
    }
    else
        error_string = "invalid instrument type";

    if (error_string != NULL)
    {
        RETURN_ERROR (error_string, "GetHeaderInput", true);
    }

#if 0
    /* Get the geo bound locations */
    this->meta.geo_bounds->min_lon = metadata->global.bounding_coords[ESPA_WEST];
    this->meta.geo_bounds->max_lon = metadata->global.bounding_coords[ESPA_EAST];
    this->meta.geo_bounds->min_lat = metadata->global.bounding_coords[ESPA_NORTH];
    this->meta.geo_bounds->max_lat = metadata->global.bounding_coords[ESPA_SOUTH];
#endif

    /* Get the map projection coordinates */
    this->meta.ul_map_corner.x = metadata->global.proj_info.ul_corner[0];
    this->meta.ul_map_corner.y = metadata->global.proj_info.ul_corner[1];
    this->meta.ul_map_corner.is_fill = true;
    this->meta.lr_map_corner.x = metadata->global.proj_info.lr_corner[0];
    this->meta.lr_map_corner.y = metadata->global.proj_info.lr_corner[1];
    this->meta.lr_map_corner.is_fill = true;

    /* Get the geo graphic coordinates */
    this->meta.ul_map_corner.x = metadata->global.ul_corner[0];
    this->meta.ul_map_corner.y = metadata->global.ul_corner[1];
    this->meta.ul_map_corner.is_fill = true;
    this->meta.lr_map_corner.x = metadata->global.lr_corner[0];
    this->meta.lr_map_corner.y = metadata->global.lr_corner[1];
    this->meta.lr_map_corner.is_fill = true;

    /* Convert the acquisition date/time values */
    sprintf (temp, "%sT%s", acq_date, acq_time);
    if (!DateInit (&this->meta.acq_date, temp, DATE_FORMAT_DATEA_TIME))
    {
        error_string = "converting acquisition date/time";
        RETURN_ERROR (error_string, "GetHeaderInput", false);
    }

    return true;
}
