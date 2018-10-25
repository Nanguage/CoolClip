from setuptools import setup, find_packages

import  re


classifiers = [
    "Development Status :: 3 - Alpha",
    "Programming Language :: Python",
    "Programming Language :: Python :: 3",
    "Programming Language :: Python :: 3.4",
    "Programming Language :: Python :: 3.5",
    "Programming Language :: Python :: 3.6",
    "Programming Language :: Python :: 3.7",
    "License :: OSI Approved :: GNU General Public License v3 (GPLv3)",
    "Intended Audience :: Science/Research",
    "Topic :: Scientific/Engineering :: Bio-Informatics",
]


keywords = [
    'genomics',
    'bioinformatics',
    'cooler',
    'Hi-C',
]


def get_version():
    with open("coolclip/__init__.py") as f:
        for line in f.readlines():
            m = re.match("__version__ = '([^']+)'", line)
            if m:
                ver = m.group(1)
                return ver
        raise IOError("Version information can not found.")


def get_long_description():
    with open("README.rst") as f:
        desc = f.read()
    return desc


def get_install_requires():
    requirements = []
    with open('requirements.txt') as f:
        for line in f:
            requirements.append(line.strip())
    return requirements


setup(
    name='coolclip',
    author='nanguage',
    author_email='nanguage@yahoo.com',
    version=get_version(),
    license='GPLv3',
    description='A small tool for clip the cooler file',
    long_description=get_long_description(),
    keywords=keywords,
    url='https://github.com/Nanguage/CoolClip',
    packages=find_packages(),
    scripts=["scripts/coolclip"],
    include_package_data=True,
    zip_safe=False,
    classifiers=classifiers,
    install_requires=get_install_requires(),
    python_requires='>=3.4, <4',
)