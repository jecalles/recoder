from setuptools import setup

with open('README.md', 'r') as handle:
    long_description = handle.read()

with open('requirements.txt') as f:
    requirements = f.read().splitlines()


setup(
    name='recoder',
    version='0.1',
    packages=[''],
    url='https://github.com/jecalles/recoder',
    license='MIT',
    author='jonathancalles',
    author_email='callesjonathan@gmail.com',
    description='a simple package for recoding overlapping genetic sequences',
    long_description=long_description,
)
