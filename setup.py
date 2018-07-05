from setuptools import setup

exec(open('main/version.py').read())

setup(
    name='ENSTWriter',
    version=__version__,
    description='',
    url='blabla',
    author='blabla',
    author_email='blabla',
    licence='MIT',
    packages=['main'],
    scripts=['bin/ENSTWriter.py'],
    zip_safe=False
)
