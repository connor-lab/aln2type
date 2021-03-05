#!/usr/bin/env python
# -*- coding: utf-8 -*-

import setuptools
from aln2type import version

requirements = [
    "pyyaml"
]

setuptools.setup(
    name="aln2type",
    version=version.__version__,
    url="https://github.com/connor-lab/aln2type",

    description="",
    long_description="",

    author="Matt Bull",
    author_email="bullmj2@gmail.com",

    maintainer="Matt Bull",
    maintainer_email="bullmj2@gmail.com",

    packages=setuptools.find_packages(),
    install_requires=requirements,

    entry_points = {
        'console_scripts': [
            'aln2type = aln2type.cli:main',
        ]
    },

    classifiers = [
        'Development Status :: 2 - Pre-Alpha',
        'Intended Audience :: Science/Research',
        'Topic :: Scientific/Engineering',
    ],

)
