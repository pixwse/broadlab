import os
import PIL
import numpy as np
# Misc helper functions


def ensure_dir_exists(path: str):
    """Make sure that the input directory exists (create it if it doesn't)
    """
    if not os.path.exists(path):
        os.makedirs(path)


def clamp(x, min_value, max_value):
    return max(min(x, max_value), min_value)


def spotname_to_xy(name: str) -> tuple[int, int]:
    # Map spot name to x,y integer coordinates
    # TODO: Error handling
    # TODO: Don't assume these spot names, should use the tissue_positions instead

    parts = name.split('_') # e.g. ['s', '002um', '00658', '01498-1']
   
    y = int(parts[2])
    parts2 = parts[3].split('-') # e.g. ['01498', '1']
    x = int(parts2[0])

    return (x,y)


def byte_barcode_to_xy(name: bytearray)  -> tuple[int, int]:
    # Convert a byte-array barcode identifier as read from the 
    # molecule_info.h5 to integer x,y positions
    
    parts = name.split(b'_')
    y = int(parts[2])
    x = int(parts[3])

    return (x, y)


def save_ndarray_as_image(input: np.ndarray, gain: float, path: str):
    scaled = input.astype(float) * gain
    scaled_u8 = scaled.clip(0, 255).astype(np.uint8)
    img_pil = PIL.Image.fromarray(scaled_u8)
    img_pil.save(path)
