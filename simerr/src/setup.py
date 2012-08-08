from distutils.core import setup
from distutils.extension import Extension
#from Cython.Distutils import build_ext

#ext_modules = [Extension("c_composite", ["c_composite.pyx"])]
ext_modules = [Extension("c_composite", ["PEassembly/c_composite.c"]),\
               Extension("c_PrimerMatch", ["PEassembly/PrimerMatch.c"]),\
		       Extension("c_tools_karkkainen_sanders", ["PEassembly/tools_karkkainen_sanders.c"])]

setup(
		name = 'PEassembly',
		version = '1.0',
		description = 'Paired-end read assembly',
		author = 'Elizabeth Tseng', 
		author_email = 'lachesis@cs.washington.edu',
		ext_package = 'PEassembly', 
		ext_modules = ext_modules,
		packages = ['PEassembly'],
		scripts = ['scripts/run_composite.py', 'scripts/composite_overlap_finder.py', 'scripts/composite_bowtie_aligner.py', 'scripts/bioconvert.py', 'scripts/remove_primers.py', 'scripts/remove_high_expected_err.py']
#		cmdclass = {'build_ext': build_ext},
 )
