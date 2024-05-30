from setuptools import setup, find_packages

setup(name='tcrconverter',
      version='0.1',
      author='Emma Bishop',
      author_email='emmab5@uw.edu',
      license='MIT',
      packages=find_packages(),
      include_package_data=True,
      package_data={'tcrconverter': ['data/*'],},
      install_requires=['pandas'])