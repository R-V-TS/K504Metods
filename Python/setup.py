from setuptools import setup, Extension

setup(
ext_modules=[Extension('Image', ['python_builder.cpp'],),Extension('DCT_filter', ['python_builder.cpp'],),Extension('DCT_frequency', ['python_builder.cpp'],)],
)
