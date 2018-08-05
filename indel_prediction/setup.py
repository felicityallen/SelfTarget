from setuptools import setup, find_packages

setup(name='predictor',
      version='0.1',
      description='Python package for Cas9 Mutation prediction',
      url='http://github.com/felicityallen/SelfTarget',
      author='Felicity Allen',
      license='MIT',
      include_package_data=True,
      packages=find_packages(),
      install_requires=['scipy','numpy>=1.9.0','matplotlib'],
      zip_safe=False)
      
      