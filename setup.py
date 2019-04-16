#!/usr/bin/env python

from setuptools import setup

def readme():
    with open('README.rst') as f:
        return f.read()


setup(
    name='qmworks',
    version='0.0.1',
    description='Automation of computations in quantum chemistry',
    license='Apache 2.0',
    url='https://github.com/sunxb05/PyFrag',
    author=["Xiaobo Sun"],
    author_email='sunxb05@gmail.com',
    keywords='chemistry workflows simulation materials',
    long_description=readme(),
    packages=["qmworks",
              "qmworks.plams",
              "qmworks.plams.core",
              "qmworks.plams.interfaces",
              "qmworks.plams.tools",
              "qmworks.components",
              "qmworks.data",
              "qmworks.data.dictionaries",
              "qmworks.hdf5",
              "qmworks.packages",
              "qmworks.parsers",
              "qmworks.templates"],
    package_data={
        "qmworks": ['data/templates/*json', 'data/dictionaries/*json', 'plams/plams_defaults']
    },
    classifiers=[
        'Intended Audience :: Science/Research',
        'Programming Language :: Python :: 3.5',
        'Development Status :: 5 - Production/Stable',
        'Intended Audience :: Science/Research',
        'Topic :: Scientific/Engineering :: Chemistry'
    ],

    install_requires=['h5py==2.6.0', 'numpy', 'noodles==0.2.3', 'pymonad',
                      'pyparsing', 'filelock', 'openpyxl', 'pyyaml', 'xlrd', 'scipy', 'tinydb', 'msgpack'
                      ]
    # pytest block
    # install_requires=['h5py', 'numpy', 'noodles==0.3.1', 'plams>=1.2', 'pymonad',
    #                   'pyparsing', 'filelock', 'openpyxl', 'pyyaml', 'xlrd', 'scipy',
    #                   'pytest', 'pytest-cov', 'pytest-mock', 'nbsphinx']
    # dependency_links=[
    #     "git+https://github.com/SCM-NV/PLAMS@master#egg=plams-1.2"]
)
