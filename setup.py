from setuptools import setup, find_packages

setup(name='tcrconvert',
      version='0.1',
      author='Emma Bishop',
      author_email='emmab5@uw.edu',
      license='MIT',
      packages=find_packages(),
      include_package_data=True,
      package_data={'tcrconvert': ['data/*'],},
      install_requires=['pandas'])