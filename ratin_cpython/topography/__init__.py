from virtualearth import VirtualEarth
try:
    from srtm import srtm
except ImportError:
    print 'SRTM loading failed'

try:
    import gdal_merge  #?gdal_merge is in the folder - no chance of ImportError
except ImportError:
    print 'GDAL not installed'