#!/usr/bin/env python

from distutils.core import setup

setup(name='kohonen',
      version='1.0',
      description='Python implementation of the Kohonen Self-Similarity Map',
      author='John Hoffman',
      author_email='jah5@princeton.edu',
      url='https://github.com/johnh2o2/kohonenpy',
      package_dir = {'kohonen': '.'},
      packages=['kohonen'],
      requires=['matplotlib', 'sklearn', 'numpy', 'scipy'],
      provides=['kohonen']
      #py_modules = [ 'vislc', 'defaults', 'vlcutils']
      #package_data={'mypkg': ['data/*.dat']},
      #packages=['pyvislc' ],
     )