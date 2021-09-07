from setuptools import setup, find_packages, Extension

setup(
    name='kmeans_CPython_API',
    version='0.0.1',
    author='Gal and Ben',
    description='CPython API',
    install_requires=['invoke',
                      'numpy',
                      'pandas',
                      'matplotlib',
                      ],
    packages=find_packages(),
    license='GPL-2',
    ext_modules=[
        # Extension(
        #     'mykmeanssp', ['kmeans.c'],
        # ),
        Extension(
            'myspkmeans', ['spkmeans.c',
                           'kmeans.c' , 'utils.c', 'goal_implementations/wam.c', 'goal_implementations/ddg.c', 'goal_implementations/lnorm.c', 'goal_implementations/spk.c'],
        )
    ],
    headers=['kmeans.h', 'spkmeans.h', 'kmeans.c', 'utils.h', 'goal_implementations/wam.h', 'goal_implementations/ddg.h', 'goal_implementations/lnorm.h', 'goal_implementations/spk.h']
)