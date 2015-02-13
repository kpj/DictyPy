from setuptools import setup


def readme():
	with open('Readme.md') as f:
		return f.read()

setup(
	name='DictyPy',
	version='0.0.3',
	description='Fun with Dictyostelium discoideum',
	long_description=readme(),
	url='https://github.com/kpj/DictyPy',
	author='kpj',
	author_email='kpjkpjkpjkpjkpjkpj@gmail.com',
	license='MIT',
	packages=[],
	test_suite='nose.collector',
	tests_require=['nose'],
	scripts=[],
	install_requires=[]
)
