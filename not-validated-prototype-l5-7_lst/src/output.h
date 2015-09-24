
#ifndef OUTPUT_H
#define OUTPUT_H


int
add_lst_band_product
(
    char *xml_filename,
    char *thermal_band_name,
    char *image_filename,
    char *product_name,
    char *band_name,
    char *short_name,
    char *long_name,
    char *data_units,
    float min_range,
    float max_range
);


#endif /* OUTPUT_H */
