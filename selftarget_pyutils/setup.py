from setuptools import setup, find_packages

setup(name='selftarget',
      version='0.1',
      description='Python utilities for Self-Target Project',
      url='http://github.com/felicityallen/SelfTarget',
      author='Felicity Allen',
      license='MIT',
      include_package_data=True,
      packages=find_packages(),
      install_requires=['scipy','numpy>=1.9.0','matplotlib'],
      zip_safe=False)
      
      