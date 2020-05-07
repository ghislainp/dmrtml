

import setuptools
from numpy.distutils.core import setup, Extension


libdmrtml = Extension(name='dmrtml_for',
                      sources=[
                            "src/dmrtml.pyf",
                            "src/czergg.f",
                            "src/dielectric_constant.f90",
                            "src/dmrtparameters.F90",
                            "src/fresnel.f90",
                            "src/soil.f90",
                            "src/disort.F90",
                            "src/dmrtml.f90",
                            "src/options.f90",
                            "src/main.f90",
                    ])


setup(
    name='DMRTML',
    version='0.1',
    author='Ghislain Picard',
    author_email='ghislain.picard@univ-grenoble-alpes.fr',
    license='GPLv3',
    classifiers=[
        "License :: OSI Approved :: MIT License",
        "Programming Language :: Python :: 3",
        "Programming Language :: Python :: 3.7",
    ],
    description='Microwave Emission Radiative Transfer Model interface',
    long_description='',
    long_description_content_type='text/markdown',
    url="https://github.com/ghislainp/dmrtml",
    ext_modules=[libdmrtml],
    packages=['dmrtml'],
    include_package_data=True,
    install_requires=['numpy'],
)
