import shutil
from distutils.core import setup
from distutils.extension import Extension
from Cython.Distutils import build_ext

#=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- core_funcs.pyx -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-#

ext_modules=[Extension("core_funcs",
             ["./src/core_funcs.pyx"],
             libraries=["m"],
             extra_compile_args = ["-ffast-math"])]

setup(
  name = "core_funcs",
  cmdclass = {"build_ext": build_ext},
  ext_modules = ext_modules)

shutil.move("./core_funcs.so", "./src/core_funcs.so")

#=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- DTD_gen.c -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=#

ext_modules=[Extension("DTD_gen", ["./src/DTD_gen.c"])]

setup(
  name = "DTD_gen",
  cmdclass = {"build_ext": build_ext},
  ext_modules = ext_modules)

shutil.move("./DTD_gen.so", "./src/DTD_gen.so")
