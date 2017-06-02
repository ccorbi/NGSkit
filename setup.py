from setuptools import setup, find_packages

setup(
    name='ngskit',
    version='0.0.1',
    author='kimlab.org',
    author_email='carles.corbi@kimlab.org',
    url = "http://kimlab.org",
    description = ("Small kit of Tools for preprocess NGS data. Customize for Kimlab's data pipelines"),
    packages=['ngskit','ngskit.utils','tests'],
    install_requires=["pandas",
                      "python-Levenshtein",
    ]



)
