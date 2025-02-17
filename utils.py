

# ----------------------------------------------------------------------------
# Helper functions

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
