from distutils.core import setup

# def readme():
#     with open('README.rst') as f:
#         return f.read()

setup(
    name='b26_toolkit',
    version='0.1.0',
    packages=['', 'src', 'src.scripts', 'src.plotting', 'src.instruments', 'src.instruments.labview_fpga_lib',
              'src.instruments.labview_fpga_lib.old', 'src.instruments.labview_fpga_lib.read_fifo',
              'src.instruments.labview_fpga_lib.galvo_scan', 'src.instruments.labview_fpga_lib.read_ai_ao',
              'src.instruments.labview_fpga_lib.pid_loop_simple',
              'src.instruments.labview_fpga_lib.labview_helper_functions', 'src.data_processing', 'tests',
              'tests.scripts', 'tests.instruments'],
    url='https://github.com/LISE-B26/b26_toolkit',
    license='GPL',
    author='Aaron Kabcenell, Jan Gieseler, and Arthur Safira',
    author_email='',
    description='Instruments, Scripts, and other classes for use with PyLabControl',
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
    tests_require=['nose'],
)
