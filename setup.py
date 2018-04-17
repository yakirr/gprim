from setuptools import setup

setup(name='gprim',
    version='0.22',
    description='Utility functions for working with statistical genetic data',
    url='http://github.com/yakirr/ypy',
    author='Yakir Reshef',
    author_email='yakir@seas.harvard.edu',
    license='MIT',
    packages=['gprim'],
    install_requires=['pysnptools','pybedtools','ypy'],
    zip_safe=False)
