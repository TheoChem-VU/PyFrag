#!/usr/bin/env python




import setuptools

with open("README.md", "r") as fh:
    long_description = fh.read()

setuptools.setup(
    name="PyFrag",
    version="0.0.1",
    author="Example Author",
    author_email="author@example.com",
    description="A small example package",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/pypa/sampleproject",
    packages=['gittest'],
    # packages=setuptools.find_packages(),
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
)




# from setuptools import setup

# def readme():
#     with open('README.rst') as f:
#         return f.read()


# setup(
#     name='qmworks',
#     version='0.3.0',
#     description='Automation of computations in quantum chemistry',
#     license='Apache 2.0',
#     url='https://github.com/sunxb05/PyFrag',
#     author=["Xiaobo Sun"]
#     author_email='f.zapata@esciencecenter.nl',
#     keywords='chemistry workflows simulation materials',
#     long_description=readme(),
#     package_dir={'': 'qmworks'},
#     packages=["qmflows",
#               "qmflows.components",
#               "qmflows.data",
#               "qmflows.data.dictionaries",
#               "qmflows.data.coskf",
#               "qmflows.examples",
#               "qmflows.examples.Conditional_workflows",
#               "qmflows.examples.Constrained_and_TS_optimizations",
#               "qmflows.examples.FDE_Fragments",
#               "qmflows.hdf5",
#               "qmflows.packages",
#               "qmflows.parsers",
#               "qmflows.templates"],
#     package_data={
#         "qmflows": ['data/templates/*json', 'data/dictionaries/*json', 'data/coskf/*coskf']
#     },
#     classifiers=[
#         'Intended Audience :: Science/Research',
#         'Programming Language :: Python :: 3.6',
#         'Development Status :: 5 - Production/Stable',
#         'Intended Audience :: Science/Research',
#         'Topic :: Scientific/Engineering :: Chemistry'
#     ],
#     install_requires=['h5py', 'numpy', 'noodles==0.3.1', 'plams>=1.2', 'pymonad',
#                       'pyparsing', 'filelock', 'openpyxl', 'pyyaml', 'xlrd', 'scipy'],
#     dependency_links=[
#         "git+https://github.com/SCM-NV/PLAMS@master#egg=plams-1.2"],

#     extras_require={
#         'test': ['pytest', 'pytest-cov', 'pytest-mock', 'nbsphinx', 'pygraphviz'],
#         'doc': ['sphinx', 'sphinx_rtd_theme', 'nbsphinx']
#     }
