from setuptools import setup, find_packages


requirements = [
            'pandas',
            'numpy',
            'scikit-learn',
            'python-Levenshtein',
            'biopython']


setup(
    name='ngskit',
    version='0.1',
    packages=  find_packages(),
    author='kimlab.org',
    author_email='carles.corbi@kimlab.org',
    url = "http://kimlab.org",
    description = ("Small kit of Tools for preprocess NGS data. Customize for Kimlab's data pipelines"),
    install_requires=requirements,
    
    classifiers=[
       "License :: OSI Approved :: MIT License",
       "Programming Language :: Python :: 3",
       "Topic :: Scientific/Engineering :: Bio-Informatics",
       ],
    license='MIT',

    entry_points={
        'console_scripts': [
                'demultiplexer = ngskit.demultiplex_reads:main',
                'trimer = ngskit.trim_reads:main',
                'barcodes = ngskit.barcodes:main',
                'counterreads = ngskit.counter_reads:main',
        ]
    }




)
