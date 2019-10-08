#!/usr/bin/env python

import os

from setuptools import (find_packages, setup)

here = os.path.abspath(os.path.dirname(__file__))

# To update the package version number, edit DeepQMC/__version__.py
version = {}
with open(os.path.join(here, 'pdb2sql', '__version__.py')) as f:
    exec(f.read(), version)

with open('README.md') as readme_file:
    readme = readme_file.read()

setup(
    name='pdb2sql',
    version=version['__version__'],
    description="PDB parser using SQL queries",
    long_description=readme + '\n\n',
    long_description_content_type='text/markdown',
    author=["Nicolas Renaud"],
    author_email='n.renaud@esciencecenter.nl',
    url='https://github.com/DeepRank/pdb2sql',
    packages=find_packages(),
    package_dir={
        'pdb2sql': 'pdb2sql'},
    include_package_data=True,
    license="Apache Software License 2.0",
    zip_safe=False,
    keywords='PDB2SQL',
    classifiers=[
            'Development Status :: 2 - Pre-Alpha',
            'Intended Audience :: Developers',
            'License :: OSI Approved :: Apache Software License',
            'Natural Language :: English',
            'Intended Audience :: Science/Research',
            'Programming Language :: Python :: 3.7',
            'Topic :: Scientific/Engineering :: Bio-Informatics'],
    test_suite='tests',
    install_requires=[
        'cython',
        'sqlalchemy',
        'matplotlib',
        'numpy',
        'schema',
        'tqdm'],
    extras_require={
        'dev': [
            'prospector[with_pyroma]',
            'yapf',
            'isort'],
        'doc': [
            'recommonmark',
            'sphinx',
            'sphinx_rtd_theme'],
        'test': [
            'coverage',
            'pycodestyle',
            'pytest',
            'pytest-cov',
            'pytest-runner'],
    })
