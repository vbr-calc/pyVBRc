#!/usr/bin/env python

"""The setup script."""

from setuptools import setup, find_packages

with open('README.md') as readme_file:
    readme = readme_file.read()

with open('requirements.txt') as reqs_file:
    requirements = [line.strip() for line in reqs_file.readlines()]

test_requirements = ['pytest>=3', ]

setup(
    author="Chris Havlin",
    author_email='chris.havlin@gmail.com',
    python_requires='>=3.10',
    classifiers=[
        'Development Status :: 2 - Pre-Alpha',
        'Intended Audience :: Developers',
        'License :: OSI Approved :: MIT License',
        'Natural Language :: English',
        'Programming Language :: Python :: 3',
        'Programming Language :: Python :: 3.10',
    ],
    description="python utilities for working with VBRc matlab output",
    entry_points={
        'console_scripts': [
            'pyVBRc=pyVBRc.cli:main',
        ],
    },
    install_requires=requirements,
    license="MIT license",
    long_description=readme,
    include_package_data=True,
    keywords='pyVBRc',
    name='pyVBRc',
    packages=find_packages(include=['pyVBRc', 'pyVBRc.*']),
    test_suite='tests',
    tests_require=test_requirements,
    url='https://github.com/chrishavlin/pyVBRc',
    version='0.1.0',
    zip_safe=False,
)
