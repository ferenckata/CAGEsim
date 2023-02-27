#!/usr/bin/env python

"""The setup script."""

from setuptools import setup, find_packages

with open('README.rst') as readme_file:
    readme = readme_file.read()

with open('HISTORY.rst') as history_file:
    history = history_file.read()

requirements = ['typer', ]

test_requirements = ['pytest>=3', "hypothesis"]

setup(
    author="Kata Ferenc",
    author_email='katalitf@uio.no',
    python_requires='>=3.6',
    classifiers=[
        'Development Status :: 2 - Pre-Alpha',
        'Intended Audience :: Developers',
        'License :: OSI Approved :: MIT License',
        'Natural Language :: English',
        'Programming Language :: Python :: 3',
        'Programming Language :: Python :: 3.6',
        'Programming Language :: Python :: 3.7',
        'Programming Language :: Python :: 3.8',
    ],
    description="Simulation code to generate CAGE-like reads with TF binding motifs",
    entry_points={
        'console_scripts': [
            'cagesim=cagesim.cli:main',
        ],
    },
    install_requires=requirements,
    license="MIT license",
    long_description=readme + '\n\n' + history,
    include_package_data=True,
    keywords='cagesim',
    name='cagesim',
    packages=find_packages(include=['cagesim', 'cagesim.*']),
    test_suite='tests',
    tests_require=test_requirements,
    url='https://github.com/ferenckata/cagesim',
    version='0.0.1',
    zip_safe=False,
)
