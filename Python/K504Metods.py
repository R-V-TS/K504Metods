# This file was automatically generated by SWIG (http://www.swig.org).
# Version 3.0.12
#
# Do not make changes to this file unless you know what you are doing--modify
# the SWIG interface file instead.

from sys import version_info as _swig_python_version_info
if _swig_python_version_info >= (2, 7, 0):
    def swig_import_helper():
        import importlib
        pkg = __name__.rpartition('.')[0]
        mname = '.'.join((pkg, '_K504Metods')).lstrip('.')
        try:
            return importlib.import_module(mname)
        except ImportError:
            return importlib.import_module('_K504Metods')
    _K504Metods = swig_import_helper()
    del swig_import_helper
elif _swig_python_version_info >= (2, 6, 0):
    def swig_import_helper():
        from os.path import dirname
        import imp
        fp = None
        try:
            fp, pathname, description = imp.find_module('_K504Metods', [dirname(__file__)])
        except ImportError:
            import _K504Metods
            return _K504Metods
        try:
            _mod = imp.load_module('_K504Metods', fp, pathname, description)
        finally:
            if fp is not None:
                fp.close()
        return _mod
    _K504Metods = swig_import_helper()
    del swig_import_helper
else:
    import _K504Metods
del _swig_python_version_info

try:
    _swig_property = property
except NameError:
    pass  # Python < 2.2 doesn't have 'property'.

try:
    import builtins as __builtin__
except ImportError:
    import __builtin__

def _swig_setattr_nondynamic(self, class_type, name, value, static=1):
    if (name == "thisown"):
        return self.this.own(value)
    if (name == "this"):
        if type(value).__name__ == 'SwigPyObject':
            self.__dict__[name] = value
            return
    method = class_type.__swig_setmethods__.get(name, None)
    if method:
        return method(self, value)
    if (not static):
        if _newclass:
            object.__setattr__(self, name, value)
        else:
            self.__dict__[name] = value
    else:
        raise AttributeError("You cannot add attributes to %s" % self)


def _swig_setattr(self, class_type, name, value):
    return _swig_setattr_nondynamic(self, class_type, name, value, 0)


def _swig_getattr(self, class_type, name):
    if (name == "thisown"):
        return self.this.own()
    method = class_type.__swig_getmethods__.get(name, None)
    if method:
        return method(self)
    raise AttributeError("'%s' object has no attribute '%s'" % (class_type.__name__, name))


def _swig_repr(self):
    try:
        strthis = "proxy of " + self.this.__repr__()
    except __builtin__.Exception:
        strthis = ""
    return "<%s.%s; %s >" % (self.__class__.__module__, self.__class__.__name__, strthis,)

try:
    _object = object
    _newclass = 1
except __builtin__.Exception:
    class _object:
        pass
    _newclass = 0


def DCT(image_, width_, height_, wind_size_):
    return _K504Metods.DCT(image_, width_, height_, wind_size_)
DCT = _K504Metods.DCT

def ADCT(image_dct, width_, height_, wind_size_):
    return _K504Metods.ADCT(image_dct, width_, height_, wind_size_)
ADCT = _K504Metods.ADCT

def getImageBlock(*args):
    return _K504Metods.getImageBlock(*args)
getImageBlock = _K504Metods.getImageBlock

def abs(zn):
    return _K504Metods.abs(zn)
abs = _K504Metods.abs

def MSE(P_image, Q_image, width, height, channels):
    return _K504Metods.MSE(P_image, Q_image, width, height, channels)
MSE = _K504Metods.MSE

def PSNR(P_image, Q_image, width, height, channels):
    return _K504Metods.PSNR(P_image, Q_image, width, height, channels)
PSNR = _K504Metods.PSNR

def PSNRHVSM(P_image, Q_image, width, height, channels):
    return _K504Metods.PSNRHVSM(P_image, Q_image, width, height, channels)
PSNRHVSM = _K504Metods.PSNRHVSM

def DCT_GPU(image, width, height, channels, window_size):
    return _K504Metods.DCT_GPU(image, width, height, channels, window_size)
DCT_GPU = _K504Metods.DCT_GPU

def ADCT_GPU(DCT_array, width, height, channels, window_size):
    return _K504Metods.ADCT_GPU(DCT_array, width, height, channels, window_size)
ADCT_GPU = _K504Metods.ADCT_GPU
TYPE_BGR = _K504Metods.TYPE_BGR
TYPE_RGB = _K504Metods.TYPE_RGB
TYPE_GRAYSCALE = _K504Metods.TYPE_GRAYSCALE
TYPE_3ARRAY = _K504Metods.TYPE_3ARRAY
class RAWImage(_object):
    __swig_setmethods__ = {}
    __setattr__ = lambda self, name, value: _swig_setattr(self, RAWImage, name, value)
    __swig_getmethods__ = {}
    __getattr__ = lambda self, name: _swig_getattr(self, RAWImage, name)
    __repr__ = _swig_repr

    def __init__(self, pixel_array, type):
        this = _K504Metods.new_RAWImage(pixel_array, type)
        try:
            self.this.append(this)
        except __builtin__.Exception:
            self.this = this
    __swig_destroy__ = _K504Metods.delete_RAWImage
    __del__ = lambda self: None

    def show(self):
        return _K504Metods.RAWImage_show(self)

    def save(self, destination_folder):
        return _K504Metods.RAWImage_save(self, destination_folder)

    def ApplyDCT(self, window_size, threshold, device):
        return _K504Metods.RAWImage_ApplyDCT(self, window_size, threshold, device)

    def DCTCoefficients(self, device):
        return _K504Metods.RAWImage_DCTCoefficients(self, device)

    def AddNoise(self, mu, sigma):
        return _K504Metods.RAWImage_AddNoise(self, mu, sigma)

    def printImageCharacteristics(self):
        return _K504Metods.RAWImage_printImageCharacteristics(self)

    def getImage(self, pixel_array):
        return _K504Metods.RAWImage_getImage(self, pixel_array)
RAWImage_swigregister = _K504Metods.RAWImage_swigregister
RAWImage_swigregister(RAWImage)

# This file is compatible with both classic and new-style classes.


