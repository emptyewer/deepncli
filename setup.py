#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""The setup script."""
import os
from setuptools import setup, find_packages

try:  # for pip >= 10
    from pip._internal import req
except ImportError:  # for pip <= 9.0.3
    from pip import req

with open('README.rst') as readme_file:
    readme = readme_file.read()

with open('HISTORY.rst') as history_file:
    history = history_file.read()

links = []
requires = []
try:
    requirements = list(req.parse_requirements(os.path.abspath('requirements.txt')))
except TypeError:
    try:
        from pip import download
    except ImportError:
        from pip._internal import download

    # new versions of pip requires a session
    requirements = req.parse_requirements(os.path.abspath('requirements.txt'), session=download.PipSession())

for item in requirements:
    # we want to handle package names and also repo urls
    if getattr(item, 'url', None):  # older pip has url
        links.append(str(item.url))
    if getattr(item, 'link', None):  # newer pip has link
        links.append(str(item.link))
    if item.req:
        requires.append(str(item.req))
print(os.listdir('.'))


setup_requirements = ['pytest-runner', ]

test_requirements = ['pytest', ]

setup(
    author="Venky Krishnamani",
    author_email='venky.krishna@icloud.com',
    classifiers=[
        'Development Status :: 2 - Pre-Alpha',
        'Intended Audience :: Developers',
        'License :: OSI Approved :: MIT License',
        'Natural Language :: English',
        'Programming Language :: Python :: 2.7',
    ],
    description="The CLI version of DEEPN bioinformatics pipeline.",
    entry_points={
        'console_scripts': [
            'deepn=deepncli.cli:main',
        ],
    },
    install_requires=requires,
    license="MIT license",
    long_description=readme + '\n\n' + history,
    include_package_data=True,
    keywords='deepncli',
    name='deepncli',
    packages=find_packages(include=['deepncli', 'deepncli.db', 'deepncli.junction', 'deepncli.utils'],
                           exclude=['deepncli.data']),
    setup_requires=setup_requirements,
    test_suite='tests',
    tests_require=test_requirements,
    url='https://github.com/emptyewer/deepncli',
    version='0.1.0',
    zip_safe=False,
)
