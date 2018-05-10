#! /usr/bin/env python

'''
    File: st_generate_qa.py

    Purpose: Builds a quality band for the surface temperature product.

    Project: Land Satellites Data Systems Science Research and Development
             (LSRD) at the USGS EROS

    License: NASA Open Source Agreement 1.3

    Algorithm author:

        Dr. Kelly G. Laraby
        Chester F. Carlson Center for Imaging Science
        Rochester Institute of Technology
'''

import os
import sys
import logging
import datetime
from argparse import ArgumentParser
from collections import namedtuple

import numpy as np
from lxml import objectify as objectify
from osgeo import gdal, osr
from scipy.ndimage.interpolation import map_coordinates

import st_utilities as util
from st_exceptions import MissingBandError
from espa import Metadata

SourceInfo = namedtuple('SourceInfo', ('proj4', 'filename'))
ThermalConstantInfo = namedtuple('ThermalInfo', ('k1', 'k2'))


def convert_radiance_to_temperature(radiance, k1, k2):
    """Calculates radiance image given temperature image using the inverse 
       Planck Function, modified to use K1 and K2 constants for a specific 
       instrument.

    Args:
        radiance <numpy.2darray>: Radiance values
        k1 <float>: K1 thermal conversion constant for the satellite
        k2 <float>: K2 thermal conversion constant for the satellite

    Returns:
        temperature <numpy.2darray>: temperature values
    """

    # Convert radiance to temperature
    return k2 / np.log(k1 / radiance + 1)


def convert_temperature_to_radiance(temperature, k1, k2):
    """Calculates temperature image given radiance image using the Planck  
       Function, modified to use K1 and K2 constants for a specific instrument.

    Args:
        temperature <numpy.2darray>: Temperature values
        k1 <float>: K1 thermal conversion constant for the satellite
        k2 <float>: K2 thermal conversion constant for the satellite

    Returns:
        radiance <numpy.2darray>: Radiance values
    """

    # Convert temperature to radiance 
    return k1 / (np.exp(k2 / temperature) - 1) 


def get_transmission_uncertainty(tau_values):
    """Calculates the transmission uncertainty term, which is part of the
       surface temperature uncertainty estimation.

    Args:
        tau_values <numpy.2darray>: Transmission values

    Returns:
        tau_uncertainty <numpy.2darray>: Transmission uncertainty values
    """

    # Just give it the right shape
    tau_uncertainty = np.empty_like(tau_values)

    # Set a lower bound, above which a quadratic fit will be made,
    # and below which last value of the fit line is extended as a constant.
    # The lower bound was set as simply the smallest transmission value from
    # the validation set that was used.
    tau_value_where_real_data_ends = 0.3005

    # Transmission coefficients calculated from MODTRAN simulations using MERRA
    tau_poly_coeff1 = -0.050527295343549
    tau_poly_coeff2 = 0.029930689405143
    tau_poly_coeff3 = 0.019127148003052

    # Calculate uncertainty values for the quadratic region
    quadratic_region = np.where(tau_values >= tau_value_where_real_data_ends)
    tau_uncertainty[quadratic_region] = \
        tau_poly_coeff1 * tau_values[quadratic_region]**2 \
        + tau_poly_coeff2 * tau_values[quadratic_region] + tau_poly_coeff3

    # Calculate uncertainty values for the constant region
    constant_region = np.where(tau_values < tau_value_where_real_data_ends)
    tau_uncertainty[constant_region] = \
        tau_poly_coeff1 * tau_value_where_real_data_ends**2 \
        + tau_poly_coeff2 * tau_value_where_real_data_ends + tau_poly_coeff3

    # Memory cleanup
    del quadratic_region
    del constant_region

    return tau_uncertainty


def get_upwelled_uncertainty(lu_values):
    """Calculates the upwelled radiance uncertainty term, which is part of the
       surface temperature uncertainty estimation.

    Args:
        lu_values <numpy.2darray>: Upwelled radiance values

    Returns:
        lu_uncertainty <numpy.2darray>: Upwelled radiance uncertainty values
    """

    # Just give it the right shape
    upwelled_uncertainty = np.empty_like(lu_values)

    # Set a lower bound, above which a quadratic fit will be made,
    # and below which last value of the fit line is extended as a constant.
    # The lower bound was set as simply the smallest upwelled radiance value
    # from the validation set that was used.
    lu_value_where_real_data_ends = 5.6788

    # Transmission coefficients calculated from MODTRAN simulations using MERRA
    lu_poly_coeff1 = -0.007147307279714
    lu_poly_coeff2 = 0.082335806862813
    lu_poly_coeff3 = -0.006782188536986

    # Calculate uncertainty values for the quadratic region
    quadratic_region = np.where(lu_values <= lu_value_where_real_data_ends)
    upwelled_uncertainty[quadratic_region] = \
        lu_poly_coeff1 * lu_values[quadratic_region]**2 \
        + lu_poly_coeff2 * lu_values[quadratic_region] + lu_poly_coeff3

    # Calculate uncertainty values for the constant region
    constant_region = np.where(lu_values > lu_value_where_real_data_ends)
    upwelled_uncertainty[constant_region] = \
        lu_poly_coeff1 * lu_value_where_real_data_ends**2 \
        + lu_poly_coeff2 * lu_value_where_real_data_ends + lu_poly_coeff3

    # Memory cleanup
    del quadratic_region
    del constant_region

    return upwelled_uncertainty


def get_downwelled_uncertainty(ld_values):
    """Calculates the downwelled radiance uncertainty term, which is part of the
       surface temperature uncertainty estimation.

    Args:
        ld_values <numpy.2darray>: Downwelled radiance values

    Returns:
        ld_uncertainty <numpy.2darray>: Downwelled radiance uncertainty values
    """

    # Just give it the right shape
    downwelled_uncertainty = np.empty_like(ld_values)

    # Set a lower bound, above which a quadratic fit will be made,
    # and below which last value of the fit line is extended as a constant.
    # The lower bound was set as simply the smallest downwelled radiance value
    # from the validation set that was used.
    ld_value_where_real_data_ends = 7.2307

    # Transmission coefficients calculated from MODTRAN simulations using MERRA
    ld_poly_coeff1 = -0.005291631300309
    ld_poly_coeff2 = 0.073763835328557
    ld_poly_coeff3 = -0.007066004902229

    # Calculate uncertainty values for the quadratic region
    quadratic_region = np.where(ld_values <= ld_value_where_real_data_ends)
    downwelled_uncertainty[quadratic_region] = \
        ld_poly_coeff1 * ld_values[quadratic_region]**2 \
        + ld_poly_coeff2 * ld_values[quadratic_region] + ld_poly_coeff3

    # Calculate uncertainty values for the constant region
    constant_region = np.where(ld_values > ld_value_where_real_data_ends)
    downwelled_uncertainty[constant_region] = \
        ld_poly_coeff1 * ld_value_where_real_data_ends**2 \
        + ld_poly_coeff2 * ld_value_where_real_data_ends + ld_poly_coeff3

    # Memory cleanup
    del quadratic_region
    del constant_region

    return downwelled_uncertainty


def get_cross_correlation(dLT_dTAU, dLT_dLU, dLT_dLD, S_TAU, S_LU, S_LD):
    """Calculates the cross correlation term, which is part of the surface
       temperature uncertainty estimation.

    Args:
        dLT_dTAU <numpy.2darray>: transmission partial
        dLT_dLU <numpy.2darray>: upwelled partial
        dLT_dLD <numpy.2darray>: downwelled partial
        S_TAU <numpy.2darray>: transmission uncertainty
        S_LU <numpy.2darray>: upwelled uncertainty
        S_LD <numpy.2darray>: downwelled uncertainty

    Returns:
        cross_correlation <numpy.2darray>: cross correlation term
    """

    # Correlation coefficients from MODTRAN simulations using MERRA.
    corr_tau_lu = -0.9899
    corr_tau_ld = -0.9857
    corr_lu_ld = 0.9965

    # Calculate cross correlation terms
    part_tau_lu = 2 * corr_tau_lu * dLT_dTAU * dLT_dLU * S_TAU * S_LU
    part_tau_ld = 2 * corr_tau_ld * dLT_dTAU * dLT_dLD * S_TAU * S_LD
    part_lu_ld = 2 * corr_lu_ld * dLT_dLU * dLT_dLD * S_LU * S_LD

    # Calculate cross correlation
    cross_correlation = part_tau_lu + part_tau_ld + part_lu_ld

    # Memory cleanup
    del part_tau_lu
    del part_tau_ld
    del part_lu_ld

    return cross_correlation


def get_unknown_uncertainty(cloud_distances, transmission_values):
    """Calculates the unknown uncertainty term, which is part of the surface
       temperature uncertainty estimation.

    Args:
        cloud_distances <numpy.2darray>: distances to nearest cloud
        transmission_values <numpy.2darray>: transmission

    Returns:
        unknown_uncertainty <numpy.2darray>: interpolated unknown uncertainty
    """

    # Flatten the inputs.
    flat_cloud_distances = cloud_distances.flatten()
    flat_transmission_values = transmission_values.flatten()

    # Matrix of "unknown errors," which was calculated from observed and
    # predicted ST errors from the L7 global validation study.
    unknown_error_matrix = np.array([[2.3749, 2.6962, 2.5620, 2.1131],
                                     [1.9912, 1.6789, 1.4471, 1.2739],
                                     [1.7925, 1.0067, 0.9143, 0.6366],
                                     [1.9416, 1.3558, 0.7604, 0.6682],
                                     [1.3861, 0.8269, 0.7404, 0.3125]])

    # tau bins are 0.3 - 0.55, 0.55 - 0.7, 0.7 - 0.85, 0.85 - 1.0
    # cloud bins are 0 - 1 km, 1 - 5 km, 5 - 10 km, 10 - 40 km, 40 - inf
    #
    # tau_interp and cloud_interp should be a vector of the center values
    # for each bin, but we also want anything outside the entire range to be
    # equal to the nearest value from the unknown matrix.
    tau_interp = np.array([0.0, 0.425, 0.625, 0.775, 0.925, 1.0])
    cld_interp = np.array([0, 0.5, 3.0, 7.5, 25.0, 82.5, 200.0])
    tau_interp_length = len(tau_interp)
    cld_interp_length = len(cld_interp)

    # Define the highest values in the vectors. These are the last values.
    tau_highest = tau_interp[-1]
    cld_highest = cld_interp[-1]

    # From input transmission values, find closest indices from tau_interp
    # vector.
    tau_close_index = np.searchsorted(tau_interp, flat_transmission_values,
                                      side='right')
    tau_close_index = tau_close_index - 1

    # Add a bin to handle values that are outside the range.
    tau_interp = np.append(tau_interp, 1.0)

    # Calculate step in vector.
    tau_step = tau_interp[tau_close_index + 1] - tau_interp[tau_close_index]

    # Calculate fractional tau index.
    with np.errstate(divide='ignore'):
        tau_frac_index = tau_close_index + ((flat_transmission_values
                                         - tau_interp[tau_close_index])
                                        / tau_step)
    high_locations = np.where(flat_transmission_values >= tau_highest)
    tau_frac_index[high_locations] = tau_interp_length - 1

    # Memory cleanup
    del tau_interp
    del tau_close_index
    del tau_step
    del high_locations

    # From input cloud distance values, find closest indices from cld_interp
    # vector.
    cld_close_index = np.searchsorted(cld_interp, flat_cloud_distances,
                                      side='right')
    cld_close_index = cld_close_index - 1

    # Add a bin to handle values that are outside the range.
    cld_interp = np.append(cld_interp, 200.0)

    # Calculate step in vector.
    cld_step = cld_interp[cld_close_index + 1] - cld_interp[cld_close_index]

    # Calculate fractional cloud index.
    with np.errstate(divide='ignore'):
        cld_frac_index = cld_close_index + ((flat_cloud_distances
                                         - cld_interp[cld_close_index])
                                         / cld_step)
    high_locations = np.where(flat_cloud_distances >= cld_highest)
    cld_frac_index[high_locations] = cld_interp_length - 1

    # Memory cleanup
    del cld_interp
    del cld_close_index
    del cld_step
    del high_locations

    # Merge the arrays so they represent coordinates in the unknown_error_matrix
    # grid to map_coordinates.
    coordinates = np.row_stack((cld_frac_index, tau_frac_index))

    # Memory cleanup
    del tau_frac_index
    del cld_frac_index

    # Interpolate the results at the specified coordinates.
    unknown_uncertainty = map_coordinates(unknown_error_matrix, coordinates,
                                          order=1, mode='nearest')

    # Memory cleanup
    del unknown_error_matrix
    del coordinates

    return unknown_uncertainty


def extract_raster_data(name, band_number):
    """Extracts raster data for the specified dataset and band number

    Args:
        name <str>: Full path dataset name
        band_number <int>: Band number for the raster to extract

    Returns:
        <raster>: 2D raster array data
    """

    dataset = gdal.Open(name)
    if dataset is None:
        raise RuntimeError('GDAL failed to open {0}'.format(name))

    raster = (dataset.GetRasterBand(band_number)
              .ReadAsArray(0, 0, dataset.RasterXSize, dataset.RasterYSize))

    return raster


def retrieve_command_line_arguments():
    """Read arguments from the command line

    Returns:
        args <arguments>: The arguments read from the command line
    """

    parser = ArgumentParser(description='Builds surface temperature QA band')

    parser.add_argument('--version',
                        action='version',
                        version=util.Version.version_text())

    parser.add_argument('--xml',
                        action='store', dest='xml_filename',
                        required=True, default=None,
                        help='The XML metadata file to use')

    parser.add_argument('--debug',
                        action='store_true', dest='debug',
                        required=False, default=False,
                        help='Output debug messages and/or keep debug data')

    args = parser.parse_args()

    # Verify that the --xml parameter was specified
    if args.xml_filename is None:
        raise Exception('--xml must be specified on the command line')

    return args


def retrieve_metadata_information(espa_metadata, band_name):
    """Reads required information from the metadata XML file

    Args:
        espa_metadata <espa.Metadata>: XML metadata
        band_name <str>: Name of band to extract metadata information from

    Returns:
        <SourceInfo>: Populated with source information
    """

    intermediate_filename = None

    # Find the intermediate band to extract information from
    for band in espa_metadata.xml_object.bands.band:
        if (band.get('product') == 'st_intermediate' and
                band.get('name') == band_name):
            intermediate_filename = str(band.file_name)

            # Get the output proj4 string
            proj4 = util.Geo.get_proj4_projection_string(intermediate_filename)

    # Error if we didn't find the required intermediate band in the data
    if intermediate_filename is None:
        raise MissingBandError('Failed to find the intermediate band'
                               ' in the input data')

    return SourceInfo(proj4=proj4, filename=intermediate_filename)


def retrieve_thermal_constants(espa_metadata, satellite):
    """Reads thermal constants from the metadata XML file

    Args:
        espa_metadata <espa.Metadata>: XML metadata
        satellite <str>: Name of satellite

    Returns:
        <ThermalConstantInfo>: Populated with thermal constants 
    """

    thermal_constants = None
    band_name = None

    # Determine the band to retrieve constants from using the satellite
    if satellite == 'LANDSAT_4' or satellite == 'LANDSAT_5':
        band_name = "b6"
    elif satellite == 'LANDSAT_7':
        band_name = "b61"
    elif satellite == 'LANDSAT_8':
        band_name = "b10"

    # Find the band to extract information from
    for band in espa_metadata.xml_object.bands.band:
        if (band.get('name') == band_name):
            thermal_constants = str(band.thermal_const)
            k1 = band.thermal_const.get('k1')
            k2 = band.thermal_const.get('k2')

    # Error if we didn't find the required intermediate band in the data
    if thermal_constants is None:
        raise MissingBandError('Failed to find the thermal band in the '
                               ' input data')

    return ThermalConstantInfo(k1=k1, k2=k2)


def calculate_qa(radiance_filename, transmission_filename, upwelled_filename,
                 downwelled_filename, emis_filename, emis_stdev_filename,
                 distance_filename, satellite, k1, k2, fill_value):
    """Calculate QA

    Args:
        radiance_filename <str>: Name of radiance file
        transmission_filename <str>: Name of atmospheric transmission file
        upwelled_filename <str>: Name of upwelled radiance file
        downwelled_filename <str>: Name of downwelled radiance file
        emis_filename <str>: Name of emissivity file
        emis_stdev_filename <str>: Name of emissivity standard deviation file
        distance_filename <str>: Name of cloud distance file
        satellite <str>: Name of satellite (e.g.: "LANDSAT_8")
        k1 <float>: K1 thermal conversion constant for the satellite
        k2 <float>: K2 thermal conversion constant for the satellite
        fill_value <float>: No data (fill) value to use

    Returns:
        <numpy.2darray>: Generated surface temperature QA band data
    """

    logger = logging.getLogger(__name__)

    logger.info('Building QA band')

    # Read the intermediate input
    # Lobs = thermal radiance
    # tau = transmission
    # Lu = upwelled radiance
    # Ld = downwelled radiance
    Lobs_array = extract_raster_data(radiance_filename, 1)
    tau_array = extract_raster_data(transmission_filename, 1)
    Lu_array = extract_raster_data(upwelled_filename, 1)
    Ld_array = extract_raster_data(downwelled_filename, 1)
    emis_array = extract_raster_data(emis_filename, 1)
    emis_stdev_array = extract_raster_data(emis_stdev_filename, 1)
    distance_array = extract_raster_data(distance_filename, 1)

    # Find fill locations.  We don't need to do this for transmission,
    # upwelled radiance, or downwelled radiance since these are created
    # with fill based on the the thermal radiance fill locations
    nonfill_locations = np.where(Lobs_array != fill_value)
    fill_locations = np.where((emis_array == fill_value) |
                              (emis_stdev_array == fill_value) |
                              (distance_array == fill_value))

    # Only operate where thermal radiance is non-fill
    Lobs = Lobs_array[nonfill_locations]
    tau = tau_array[nonfill_locations]
    Lu = Lu_array[nonfill_locations]
    Ld = Ld_array[nonfill_locations]
    emis = emis_array[nonfill_locations]
    emis_stdev = emis_stdev_array[nonfill_locations]
    distance = distance_array[nonfill_locations]

    # Memory cleanup
    del tau_array
    del Lu_array
    del Ld_array
    del emis_array
    del emis_stdev_array
    del distance_array

    # Calculate partials
    dLT_dTAU = (Lu - Lobs) / (emis * tau**2)
    dLT_dLU = -1 / (tau * emis)
    dLT_dLD = (emis - 1) / emis
    dLT_dLOBS = 1 / (tau * emis)

    # Create radiance image that will be used as nominal values around which
    # the uncertainty values will operate
    Le_Image = (Lobs - Lu - (1 - emis) * Ld * tau) / (emis * tau)
    Le_Image[Le_Image > 30.0] = 0

    # Memory cleanup
    del Lobs
    del emis

    # Calculate transmission, upwelled radiance, and downwelled radiance
    # uncertainty
    S_TAU = get_transmission_uncertainty(tau)
    S_LU = get_upwelled_uncertainty(Lu)
    S_LD = get_downwelled_uncertainty(Ld)

    # Memory cleanup
    del Lu
    del Ld

    # Calculate atmosphere uncertainty
    S_A = (dLT_dTAU * S_TAU)**2 + (dLT_dLU * S_LU)**2 + (dLT_dLD * S_LD)**2

    # Look up satellite thermal band uncertainty in radiance based on satellite.
    # These are based on the following uncertainty values in K:
    # Landsat 4: 0.8
    # Landsat 5: 0.8
    # Landsat 7: 0.4
    # Landsat 8: 0.3
    if satellite == 'LANDSAT_4':
        landsat_uncertainty = 0.08959421
    elif satellite == 'LANDSAT_5':
        landsat_uncertainty = 0.08959421
    elif satellite == 'LANDSAT_7':
        landsat_uncertainty = 0.04471454
    elif satellite == 'LANDSAT_8':
        landsat_uncertainty = 0.03577072
    else:
        raise Exception('Unsupported satellite sensor')

    # Calculate instrument uncertainty
    S_I = (dLT_dLOBS * landsat_uncertainty)**2

    # Memory cleanup
    del dLT_dLOBS

    # Get the RMSE of the linear regression fit of the spectral emissivity
    # adjustment procedure (which is the emis_data calculation in
    # estimate_landsat_emissivity)
    if satellite == 'LANDSAT_4':
        emis_regfit = 0.00085135
    elif satellite == 'LANDSAT_5':
        emis_regfit = 0.0013
    elif satellite == 'LANDSAT_7':
        emis_regfit = 0.0011
    elif satellite == 'LANDSAT_8':
        emis_regfit = 0.00093909
    else:
        raise Exception('Unsupported satellite sensor')

    # Calculate the total algorithmic uncertainty of the Temperature Emissivity
    # Separation algorithm used to produce the ASTER GED emissivities described
    # in Hulley et al. 2012
    eret13 = 0.0164
    eret14 = 0.0174
    eret = np.sqrt((eret13**2 + eret14**2) / 2)

    # Calculate emissivity uncertainty
    S_E = np.sqrt((emis_stdev**2 + emis_regfit**2 + eret**2) / 3)

    # Memory cleanup
    del emis_stdev

    # Calculate cross correlation terms
    S_P = get_cross_correlation(dLT_dTAU, dLT_dLU, dLT_dLD, S_TAU, S_LU, S_LD)

    # Memory cleanup
    del S_TAU
    del S_LU
    del S_LD
    del dLT_dTAU
    del dLT_dLU
    del dLT_dLD

    # Calculate unknown uncertainty in K
    unknown_uncertainty = get_unknown_uncertainty(distance, tau)
    S_U = unknown_uncertainty**2

    # Memory cleanup
    del tau
    del distance
    del unknown_uncertainty

    # The conversion to residual radiance must be about a nominal temperature,
    # the LST image
    Delta_Temps = S_U / 2.0

    # Memory Cleanup
    del S_U

    # Create temperature image that will be used as nominal values around which
    # the uncertainty values will operate
    LST_Image = convert_radiance_to_temperature(Le_Image, k1, k2) 

    # Find the deltas around the nominal values
    LST_Minus_Delta = LST_Image - Delta_Temps
    LST_Plus_Delta = LST_Image + Delta_Temps

    # Memory Cleanup
    del LST_Image
    del Delta_Temps

    # Convert from temperature (K) to radiance
    S_U_Radiance_High = convert_temperature_to_radiance(LST_Plus_Delta, k1, k2)
    S_U_Radiance_Low = convert_temperature_to_radiance(LST_Minus_Delta, k1, k2)

    # Memory Cleanup
    del LST_Minus_Delta
    del LST_Plus_Delta

    # The total uncertainty is the difference between the high and low values
    S_U_Radiance = S_U_Radiance_High - S_U_Radiance_Low
    S_U_Radiance[S_U_Radiance < 0.0] = 0

    # Memory Cleanup
    del S_U_Radiance_High 
    del S_U_Radiance_Low 

    # Calculate uncertainty in radiance units
    st_uncertainty_radiance = np.sqrt(S_A + S_I + S_E + S_P + S_U_Radiance)

    # Memory cleanup
    del S_A
    del S_I
    del S_E
    del S_P
    del S_U_Radiance

    # Again the st_uncertainty values are residuals about some nominal radiance
    # values, so we must convert to residual temperature about some nominal 
    # reference radiance, which we get from the original radiance image
    Delta_Radiance = st_uncertainty_radiance / 2.0

    # Find the deltas around the nominal values
    Radiance_Minus_Delta = Le_Image - Delta_Radiance
    Radiance_Plus_Delta = Le_Image + Delta_Radiance

    # Memory Cleanup
    del Le_Image 
    del Delta_Radiance 

    # Convert from radiance to temperature (K)
    Temp_uncertainty_High = convert_radiance_to_temperature(Radiance_Plus_Delta,
                                                            k1, k2)
    Temp_uncertainty_Low = convert_radiance_to_temperature(Radiance_Minus_Delta,
                                                           k1, k2)

    # Memory Cleanup
    del Radiance_Plus_Delta 
    del Radiance_Minus_Delta 

    # The total uncertainty is the difference between the high and low values
    Temp_Uncertainty = Temp_uncertainty_High - Temp_uncertainty_Low
    Temp_Uncertainty[Temp_Uncertainty > 100.0] = 0

    # Memory Cleanup
    del Temp_uncertainty_High 
    del Temp_uncertainty_Low 

    # Give st_uncertainty the same dimensions as the original Lobs
    st_uncertainty_array = Lobs_array.copy()
    st_uncertainty_array.fill(fill_value)
    st_uncertainty_array[nonfill_locations] = Temp_Uncertainty

    # Memory cleanup
    del nonfill_locations
    del Lobs_array
    del Temp_Uncertainty

    # Apply fill to results where other inputs were fill
    st_uncertainty_array[fill_locations] = fill_value

    # Memory cleanup
    del fill_locations

    return st_uncertainty_array


def add_qa_band_to_xml(espa_metadata, filename, sensor_code, no_data_value):
    """Adds the QA band to the Metadata XML file

    Args:
        espa_metadata <espa.Metadata>: XML metadata information
        filename <str>: Full path for the output file to create
        sensor_code <str>: Name prefix for the sensor
        no_data_value <float>: Value to use for fill
    """

    logger = logging.getLogger(__name__)

    logger.info('Adding {0} to Metadata XML'.format(filename))

    # Create an element maker
    maker = objectify.ElementMaker(annotate=False, namespace=None, nsmap=None)

    source_product = 'toa_refl'

    # Find the TOA Band 1 to use for the specific band details
    base_band = None
    for band in espa_metadata.xml_object.bands.band:
        if (band.get('product') == source_product and
                band.get('name') == 'toa_band1'):
            base_band = band

    if base_band is None:
        raise MissingBandError('Failed to find the TOA band 1'
                               ' in the input data')

    qa_band = maker.band()
    qa_band.set('product', 'st_qa')
    qa_band.set('source', source_product)
    qa_band.set('name', 'st_qa')
    qa_band.set('category', 'qa')
    qa_band.set('data_type', 'INT16')
    qa_band.set('scale_factor', str(SCALE_FACTOR))
    qa_band.set('nlines', base_band.attrib['nlines'])
    qa_band.set('nsamps', base_band.attrib['nsamps'])
    qa_band.set('fill_value', str(no_data_value))

    qa_band.short_name = maker.element('{0}STQA'.format(sensor_code))

    qa_band.long_name = maker.element('Surface temperature quality band')
    qa_band.file_name = maker.element(filename)

    qa_band.pixel_size = base_band.pixel_size

    qa_band.resample_method = maker.element('none')
    qa_band.data_units = maker.element('temperature (kelvin)')

    qa_band.valid_range = maker.element()
    qa_band.valid_range.set('min', '0')
    qa_band.valid_range.set('max', '32767')

    qa_band.app_version = maker.element(util.Version.app_version())

    # Get the production date and time in string format
    # Strip the microseconds and add a Z
    date_now = ('{0}Z'.format(datetime.datetime.now()
                              .strftime('%Y-%m-%dT%H:%M:%S')))
    qa_band.production_date = maker.element(date_now)

    # Append the band to the XML
    espa_metadata.xml_object.bands.append(qa_band)

    # Validate the XML
    espa_metadata.validate()

    # Write it to the XML file
    espa_metadata.write()


def write_qa_product(samps, lines, transform, wkt, no_data_value, filename,
                     file_data):
    """Creates the QA band file

    Args:
        samps <int>: Samples in the data
        lines <int>: Lines in the data
        transform <2x3:float>: GDAL Affine transformation matrix
                               [0] - Map X of upper left corner
                               [1] - Pixel size in X direction
                               [2] - Y rotation
                               [3] - Map Y of upper left corner
                               [4] - X rotation
                               [5] - Pixel size in Y direction
        wkt <str>: Well-Known-Text describing the projection
        no_data_value <float>: Value to use for fill
        filename <str>: Full path for the output file to create
        file_data <numpy.2darray>: Binary image data to be output
    """

    # Scale the data
    file_data[file_data != no_data_value] *= MULT_FACTOR

    logger = logging.getLogger(__name__)

    logger.info('Creating {0}'.format(filename))
    util.Geo.generate_raster_file(gdal.GetDriverByName('ENVI'),
                                  filename,
                                  file_data,
                                  samps,
                                  lines,
                                  transform,
                                  wkt,
                                  no_data_value,
                                  gdal.GDT_Int16)

    hdr_filename = filename.replace('.img', '.hdr')
    logger.info('Updating {0}'.format(hdr_filename))
    util.Geo.update_envi_header(hdr_filename, no_data_value)

    # Remove the *.aux.xml file generated by GDAL
    aux_filename = filename.replace('.img', '.img.aux.xml')
    if os.path.exists(aux_filename):
        os.unlink(aux_filename)


def generate_qa(xml_filename, no_data_value):
    """Provides the main processing algorithm for generating the QA product.

    Args:
        xml_filename <str>: Filename for the ESPA Metadata XML
        no_data_value <float>: No data (fill) value to use
    """

    logger = logging.getLogger(__name__)

    # XML metadata
    espa_metadata = Metadata(xml_filename)
    espa_metadata.parse()

    radiance_src_info \
        = retrieve_metadata_information(espa_metadata,
                                        'st_thermal_radiance')
    transmission_src_info \
        = retrieve_metadata_information(espa_metadata,
                                        'st_atmospheric_transmittance')
    upwelled_src_info \
        = retrieve_metadata_information(espa_metadata,
                                        'st_upwelled_radiance')
    downwelled_src_info \
        = retrieve_metadata_information(espa_metadata,
                                        'st_downwelled_radiance')
    emis_src_info = retrieve_metadata_information(espa_metadata, 'emis')
    emis_stdev_src_info \
        = retrieve_metadata_information(espa_metadata, 'emis_stdev')
    satellite = espa_metadata.xml_object.global_metadata.satellite
    thermal_info = retrieve_thermal_constants(espa_metadata, satellite)

    # Determine output information.  Make it like the emissivity band
    sensor_code = get_satellite_sensor_code(xml_filename)
    dataset = gdal.Open(emis_src_info.filename)
    output_srs = osr.SpatialReference()
    output_srs.ImportFromWkt(dataset.GetProjection())
    output_transform = dataset.GetGeoTransform()
    samps = dataset.RasterXSize
    lines = dataset.RasterYSize
    del dataset

    # Build cloud distance filename
    distance_img_filename = ''.join([xml_filename.split('.xml')[0],
                                     '_st_cloud_distance', '.img'])

    # Build QA information in memory
    qa = calculate_qa(radiance_src_info.filename,
                      transmission_src_info.filename,
                      upwelled_src_info.filename,
                      downwelled_src_info.filename,
                      emis_src_info.filename,
                      emis_stdev_src_info.filename,
                      distance_img_filename,
                      satellite,
                      float(thermal_info.k1),
                      float(thermal_info.k2),
                      no_data_value)

    # Build QA filename
    qa_img_filename = ''.join([xml_filename.split('.xml')[0],
                               '_st_uncertainty', '.img'])

    # Write QA product
    write_qa_product(samps=samps,
                     lines=lines,
                     transform=output_transform,
                     wkt=output_srs.ExportToWkt(),
                     no_data_value=no_data_value,
                     filename=qa_img_filename,
                     file_data=qa)

    add_qa_band_to_xml(espa_metadata=espa_metadata,
                       filename=qa_img_filename,
                       sensor_code=sensor_code,
                       no_data_value=no_data_value)


def get_satellite_sensor_code(xml_filename):
    """Derives the satellite-sensor code from the XML filename

    Args:
        xml_filename <str>: Filename for the ESPA Metadata XML

    Returns:
        <str>: Satellite sensor code
    """

    collection_prefixes = ['LT04', 'LT05', 'LE07', 'LT08', 'LC08', 'LO08']

    base_name = os.path.basename(xml_filename)

    satellite_sensor_code = base_name[0:4]
    if satellite_sensor_code in collection_prefixes:
        return satellite_sensor_code

    raise Exception('Satellite-Sensor code ({0}) not understood'
                    .format(satellite_sensor_code))



# Specify the no data, scale factor, and multiplication factor.
NO_DATA_VALUE = -9999
SCALE_FACTOR = 0.01
MULT_FACTOR = 100.0


def main():
    """Main processing for building the surface temperature QA band
    """

    # Command Line Arguments
    args = retrieve_command_line_arguments()

    # Check logging level
    logging_level = logging.INFO
    if args.debug:
        logging_level = logging.DEBUG

    # Setup the default logger format and level.  Log to STDOUT.
    logging.basicConfig(format=('%(asctime)s.%(msecs)03d %(process)d'
                                ' %(levelname)-8s'
                                ' %(filename)s:%(lineno)d:'
                                '%(funcName)s -- %(message)s'),
                        datefmt='%Y-%m-%d %H:%M:%S',
                        level=logging_level,
                        stream=sys.stdout)
    logger = logging.getLogger(__name__)

    logger.info('*** Begin ST Generate QA ***')

    try:
        # Register all the gdal drivers
        gdal.AllRegister()

        # Call the main processing routine
        generate_qa(xml_filename=args.xml_filename,
                    no_data_value=NO_DATA_VALUE)

    except Exception:
        logger.exception('Processing failed')
        sys.exit(1)  # EXIT FAILURE

    logger.info('*** ST Generate QA - Complete ***')

if __name__ == '__main__':
    main()
