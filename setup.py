from distutils.core import setup

setup(
    name='ko_contributions',
    version='0.1.dev',
    packages=['ko_contrib' , 'ko_contrib/scripts' , 'ko_contrib/src'],
    install_requires = ['biom-format' , 'pandas' , 'seaborn'],
    package_dir={'': 'ko_contrib'},
    url='',
    license='',
    author='xaradrim',
    author_email='',
    description='version test'
)
