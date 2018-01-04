from setuptools import setup

setup(name='mutfreq',
      version='0.1',
      description='experimental tool to analyze mutation data',
      url='https://github.com/cfblaeb/mutfreq',
      author='Lasse Ebdrup Pedersen',
      author_email='https://orcid.org/0000-0002-6064-919X',
      license='none',
      packages=['mutfreq'],
      install_requires=[
          'pandas',
      ],
      zip_safe=False)
