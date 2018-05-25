from distutils.core import setup

# def readme():
#     with open('README.md') as f:
#         return f.read()

setup(
    name='b26_toolkit',
    version='0.1.0',
    package_dir={'b26_toolkit': ''},
    packages=['b26_toolkit.pylabcontrol', 'b26_toolkit.pylabcontrol.data_processing', 'b26_toolkit.pylabcontrol.instruments', 'b26_toolkit.pylabcontrol.plotting', 'b26_toolkit.pylabcontrol.scripts', 'b26_toolkit.tests'],
    url='https://github.com/LISE-B26/b26_toolkit',
    license='GPL',
    author='Aaron Kabcenell, Jan Gieseler, and Arthur Safira',
    author_email='',
    description='Instruments, Scripts, and other classes for use with pylabcontrol',
    # long_description=readme(),
    keywords='laboratory control',
    classifiers=[
        'Programming Language :: Python :: 2.7',
        'License :: OSI Approved :: GNU General Public License v3 (GPLv3)',
        'Development Status :: 4 - Beta',
        'Environment :: Win32 (MS Windows)',
        'Intended Audience :: Science/Research',
        'Topic :: Scientific/Engineering :: Physics',
    ],
    install_requires=[
        'pillow'
    ],
    test_suite='nose.collector',
    tests_require=['nose']
)
