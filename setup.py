from distutils.core import setup
from distutils.extension import Extension
from Cython.Distutils import build_ext

ext_modules = [Extension("c_PrimerMatch", ["PrimerMatch.pyx"]),\
		Extension("c_tools_karkkainen_sanders", ["tools_karkkainen_sanders.pyx"])]

setup(
		name = 'Cython version of composite',
		cmdclass = {'build_ext': build_ext},
		ext_modules = ext_modules
 )
