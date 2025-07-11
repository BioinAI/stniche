from setuptools import setup, find_packages

setup(
    name='stNiche',
    version='0.1.0',
    packages=find_packages(),
    install_requires=[
        'pandas>=2.0.0,<2.1',
        'numpy>=1.22,<1.23',
        'scipy>=1.10,<1.11',
        'statsmodels>=0.14,<0.15',
        'matplotlib>=3.7,<3.8',
        'seaborn>=0.13,<0.14',
        'tqdm>=4.67,<4.68',
        'scanpy>=1.9.8,<1.10',
        'squidpy>=1.2,<1.3',
        'plotly>=6.0,<6.1',
        'gseapy>=1.0,<1.1',
        'kaleido>=0.2.1,<0.3'
       
    ],
    author='Mintian Cui',
    author_email='1308318910@qq.com',
    description='Spatial structure analysis toolkit for transcriptomics.',
    long_description=open('README.md', encoding='utf-8').read(),
    long_description_content_type='text/markdown',
    license='MIT',
    url='https://github.com/BioinAI/stniche',
    python_requires='>=3.9',
    classifiers=[
        'Programming Language :: Python :: 3',
        'License :: OSI Approved :: MIT License',
        'Operating System :: OS Independent',
        'Topic :: Scientific/Engineering :: Bio-Informatics',
    ],
)