
from setuptools import setup, find_packages

from codecs import open
from os import path

here=path.abspath(path.dirname(__file__))

with open(path.join(here, 'README.md'), encoding='utf-8') as f:
    long_description=f.read()

setup(
    name='extractor_hyperspectral',
    version='1.0.0',

    description='This extractor processes HDF files into netCDF.',
    long_description=long_description,

    entry_points={
        'console_scripts': [
            'terra_hyperspectral=hyperspectral.hyperspectral:main',
        ],
    },

    scripts=['bin/hyperspectral_workflow.sh',
             'bin/hyperspectral_calibration_reduction.sh',
            ],

    install_requires=[
        'enum34',
        'pyyaml',
        'pika>=0.10.0',
        'urllib3',
        'requests>=2.11.0',
        'numpy',
        'netCDF4',
        'pyclowder>=2.0.0',
        ],

    dependency_links=['https://opensource.ncsa.illinois.edu/bitbucket/rest/archive/latest/projects/CATS/repos/pyclowder2/archive?format=zip#egg=pyclowder-2.0.0'],

    packages=find_packages(),

    # basic package metadata
    url='https://github.com/terraref/extractors-hyperspectral',
    author='Max Burnette',
    author_email='mburnet2@illinois.edu',

    license='NCSA',
    classifiers=[
        'Development Status :: 3 - Alpha',
        'License :: OSI Approved :: NCSA License',

        'Programming Language :: Python :: 2',
        'Programming Language :: Python :: 2.7',
    ],
    keywords='terraref clowder extractor'

)

