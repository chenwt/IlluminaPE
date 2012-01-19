from distutils.core import setup
from distutils.extension import Extension
from Cython.Distutils import build_ext

ext_modules = [Extension("c_composite", ["c_composite.pyx"])]

setup(
		name = 'Cython version of composite',
		cmdclass = {'build_ext': build_ext},
		ext_modules = ext_modules
 )
