from setuptools import setup

setup(
    name='esqg',
    version='1.0',
    packages=['esqg'],
    entry_points={
        'console_scripts': [
            'esqg = esqg.esqg_data:main_function',
        ],
    },
    author='Tianshi Du',
    author_email='dutianshi@hotmail.com',
    description='effective surface quasi-geostrophic theory',
    url='https://github.com/DuTianshi/esqg/tree/main',
    classifiers=[
        'Programming Language :: Python :: 3',
    ],
)
