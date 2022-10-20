#!/usr/bin/env python

"""The setup script."""

from setuptools import find_packages, setup

with open("README.md") as readme_file:
    readme = readme_file.read()

with open("requirements.txt") as reqs_file:
    requirements = [line.strip() for line in reqs_file.readlines()]

test_requirements = [
    "pytest>=3",
]

setup(
    author="Chris Havlin",
    author_email="chris.havlin@gmail.com",
    python_requires=">=3.9",
    classifiers=[
        "Development Status :: 2 - Pre-Alpha",
        "Intended Audience :: Developers",
        "License :: OSI Approved :: MIT License",
        "Natural Language :: English",
        "Programming Language :: Python :: 3",
        "Programming Language :: Python :: 3.9",
        "Programming Language :: Python :: 3.10",
    ],
    description="python utilities for working with the VBRc",
    entry_points={
        "console_scripts": [
            "pyVBRc=pyVBRc.cli:main",
        ],
    },
    install_requires=requirements,
    license="MIT license",
    long_description="python utilities for working with the Very Broadband "
    "Rheology Calculator, see https://github.com/vbr-project/pyVBRc",
    include_package_data=True,
    package_data={"": ["sample_data/*.mat"]},
    keywords="pyVBRc",
    name="pyVBRc",
    packages=find_packages(include=["pyVBRc", "pyVBRc.*"]),
    test_suite="tests",
    tests_require=test_requirements,
    url="https://github.com/vbr-project/pyVBRc",
    version="0.1.1",
    zip_safe=False,
)
