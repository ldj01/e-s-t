


FUNCTION MAKE_L7_ZEROS_MASK, band6 

  image_size = SIZE(band6)
  find_zeros = WHERE(band6 EQ 0)

  zeros_mask = MAKE_ARRAY(image_size[1],image_size[2])
  mask1 = REFORM(zeros_mask[*,*])
  mask1[*] = 1 
  mask1[find_zeros] = 0

  zeros_mask[*,*] = mask1

RETURN, zeros_mask

END
