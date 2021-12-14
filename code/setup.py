# -*- coding: utf-8 -*-

from setuptools import setup, find_packages

setup(
    name='seGMM',
    version='1.2.1',
    description=(
        'seGMM is a tool that determines the gender of a sample from the called genotype data integrated BAM files and jointly considers information on the X and Y chromosomes in diverse genotype data including panel data. seGMM apply Gaussian Mixture Model (GMM) clustering to classify the samples into two clusters.'
    ),
    long_description=open('README.md').read(),
    long_description_content_type="text/markdown",
    author='Sihan Liu',
    author_email='liusihan@wchscu.cn',
    maintainer='Sihan Liu',
    maintainer_email='liusihan@wchscu.cn',
    license='MIT License',
    platforms=["linux"],
    url='https://github.com/liusihan/seGMM',
    packages=find_packages(),
    entry_points={
        'console_scripts': [
            'seGMM = seGMM.main:main',
        ]
    },
    package_data={
        'seGMM': ['script/*','data/*']
    },
    classifiers=[
        'Operating System :: OS Independent',
        'License :: OSI Approved :: MIT License',
        'Programming Language :: Python :: 3'
    ],
    python_requires='>=3'
)
