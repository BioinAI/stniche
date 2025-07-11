from setuptools import setup, find_packages

setup(
    name='stniche',
    version='0.1.0',
    packages=find_packages(),
    install_requires=[
        'pandas',
        'numpy',
        'scipy',
        'statsmodels',
        'anndata',
        'tqdm'
    ],
    author='Mintian Cui',
    author_email='1308318910@qq.com',
    description='Spatial structure analysis toolkit for transcriptomics.',
    license='MIT',
    url='https://github.com/BioinAI/stniche',
)
