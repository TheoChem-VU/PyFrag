#!/usr/bin/env python




# import setuptools

# with open("README.md", "r") as fh:
#     long_description = fh.read()

# setuptools.setup(
#     name="PyFrag",
#     version="0.0.1",
#     author="Example Author",
#     author_email="author@example.com",
#     description="A small example package",
#     long_description=long_description,
#     long_description_content_type="text/markdown",
#     url="https://github.com/pypa/sampleproject",
#     packages=['gittest'],
#     # packages=setuptools.find_packages(),
#     classifiers=[
#         "Programming Language :: Python :: 3",
#         "License :: OSI Approved :: MIT License",
#         "Operating System :: OS Independent",
#     ],
# )




from setuptools import setup

def readme():
    with open('README.md') as f:
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
    package_dir={'': 'qmworks'},
    packages=["qmworks",
              "qmworks.components",
              "qmworks.data",
              "qmworks.data.dictionaries",
              "qmworks.hdf5",
              "qmworks.packages",
              "qmworks.parsers",
              "qmworks.templates"],
    package_data={
        "qmworks": ['data/templates/*json', 'data/dictionaries/*json']
    },
    classifiers=[
        'Intended Audience :: Science/Research',
        'Programming Language :: Python :: 3.6',
        'Development Status :: 5 - Production/Stable',
        'Intended Audience :: Science/Research',
        'Topic :: Scientific/Engineering :: Chemistry'
    ]
)
