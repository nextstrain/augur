from collections import namedtuple

__all__ = ['__version__', 'version_info']

__vinfo = namedtuple('version_info', ['major', 'minor', 'micro'])

__version__ = '6.4.3'
version_info = __vinfo(*map(int, __version__.split('.')))

del __vinfo
