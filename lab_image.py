import os
import matplotlib.pyplot as plt

import PIL
import cv2
import large_image
import numpy as np
import io

# Lab stuff to work with the H&E image only - checking how we can load it, 
# split it into more convenient tiles etc.


def lab_load_large_image():
    path = '/media/erik/T9/HD_images/run2_A/2024-08-29SRParryRichterPD1460124.tif'

    # Doesn't work using directly, says:
    # Image size (9954918400 pixels) exceeds limit of 178956970 pixels
    # (meaning that it's around 100k x 100k large)
    # img = PIL.Image.open(path)

    # Testing with OpenCV
    # Only supports up to 2^30, and can't be configure to more than that
    # img = cv2.imread(path)

    # Testing large-image (with pip install large-image[tiff])
    ts = large_image.getTileSource(path)
    qq = ts.getMetadata()

    plt.interactive(True)

    for tile_info in ts.tileIterator(
            region=dict(left=40000, top=40000, width=20000, height=20000, units='base_pixels'),
            scale=dict(magnification=20),
            tile_size=dict(width=1000, height=1000),
            tile_overlap=dict(x=50, y=50),
            format=large_image.tilesource.TILE_FORMAT_PIL
            ):

        im_tile = np.array(tile_info['tile'])        
        plt.imshow(im_tile)
        plt.show()
        print('One tile ok')
        plt.waitforbuttonpress()


def split_large_image():
    # For example: tile 166 is within the ST region (top part of the
    # large crack-like structure)

    path = '/media/erik/T9/HD_images/run2_A/2024-08-29SRParryRichterPD1460124.tif'
    output_path = '_temp/cll2a_tiles_png'

    # Testing large-image (with pip install large-image[tiff])
    ts = large_image.getTileSource(path)
    qq = ts.getMetadata()

    if False:
        region, mime_type = ts.getRegion(
            output={'maxWidth': 10000, 'maxHeight': 10000, 'format': 'numpy'}
        )

        image = PIL.Image.open(io.BytesIO(region))
        im.save(os.path.join(output_path, 'lowres.jpg'))

    tile_no = 0
    for tile_info in ts.tileIterator(
            # region=dict(left=0, top=0, width=20000, height=20000, units='base_pixels'),
            # scale=dict(magnification=20),
            tile_size=dict(width=5000, height=5000),
            tile_overlap=dict(x=0, y=0),
            format=large_image.tilesource.TILE_FORMAT_PIL
            ):

        if tile_no > 160 and tile_no < 170:
            x = tile_info['x']
            y = tile_info['y']
            im_tile = np.array(tile_info['tile'])        
            im = PIL.Image.fromarray(im_tile)
            im.save(os.path.join(output_path, f'tile_{x}_{y}.png'))
            print(f'Tile {tile_no} saved')

        tile_no += 1


if __name__ == '__main__':
    # lab_load_large_image()
    split_large_image()
