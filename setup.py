#!/usr/bin/env python3
# -*- encoding: utf-8 -*-

import os
from setuptools import setup
import confilter


def read(*names):
    values = {}
    extensions = ['', '.txt', '.rst', '.md',]
    for name in names:
        v = ''
        for ext in extensions:
            filename = name + ext
            if os.path.isfile(filename):
                with open(filename) as f:
                    v = f.read()
        values[name] = v
    return values


long_description = """%(README)s""" % read('README')


setup(
    name='confilter',
    version=confilter.VERSION,
    description='Conformer Filtration',
    long_description=long_description,
    long_description_content_type='text/markdown',
    classifiers=[
        "Development Status :: 3 - Alpha",
        "Environment :: Console",
        "Intended Audience :: Developers",
        "Intended Audience :: Science/Research",
        "Programming Language :: Python :: 3",
        "Topic :: Software Development :: Build Tools",
        "Topic :: Scientific/Engineering :: Chemistry",
        "Topic :: Scientific/Engineering :: Information Analysis",
        "Topic :: Scientific/Engineering :: Visualization",
    ],
    author='Orlando Acevedo',
    author_email='orlando.acevedo@miami.edu',
    keywords='Conformer Filtration',
    maintainer='Xiang Zhong',
    maintainer_email='xxz385@miami.edu',
    url='https://github.com/orlandoacevedo/confilter',
    license='MIT',
    packages=['confilter'],
    entry_points={
        'console_scripts': [
            'confilter = confilter.main:parsecmd',
        ]
    },
    install_requires=['matplotlib'],
)



