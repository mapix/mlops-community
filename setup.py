from setuptools import setup, find_packages

install_requires = []

setup(
    name='mlops_helloworld',
    version='0.0.1',
    author='dp.tech',
    author_email='mlops@dp.tech',
    description=('Project-level MLOps'),
    url='',
    license=None,
    keywords='MLOps Hello World',
    install_requires=install_requires,
    packages=find_packages(),
    zip_safe=False,
    # packages=packages,
    entry_points={'console_scripts': [
        'benchmark_mlops = mlops_helloworld.benchmark:main',
    ]},
    include_package_data=True)


setup(
    name='benchmark_helloworld',
    version='0.0.1',
    author='dp.tech',
    author_email='mlops@dp.tech',
    description=('Benchmark helloworld.'),
    url='',
    license=None,
    keywords='MLOps Hello World',
    install_requires=install_requires,
    packages=find_packages(),
    zip_safe=False,
    # packages=packages,
    entry_points={'console_scripts': [
        'benchmark_demo_ = benchmark_helloworld.benchmark:main',
    ]},
    include_package_data=True)
