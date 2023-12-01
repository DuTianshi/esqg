from setuptools import setup

setup(
    name='esqg',
    version='1.0',
    py_modules=['esqg'],
    install_requires=[
        'numpy',
        'seawater',
        'scikit-learn',
        # Add other dependencies here
    ],
    entry_points={
        'console_scripts': [
            'esqg-script = esqg:main',
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
