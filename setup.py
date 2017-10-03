from setuptools import setup

setup(name='mutfreq',
      version='0.1',
      description='sam to df',
      url='https://github.com/cfblaeb/mutfreq',
      author='Lasse Ebdrup Pedersen',
      author_email='gotodaplace@hotmail.com',
      license='none',
      packages=['mutfreq'],
      install_requires=[
          'pandas',
      ],
      zip_safe=False)
