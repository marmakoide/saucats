#! /usr/bin/env python

# the following two variables are used by the target "waf dist"
VERSION = '1.0.0'
APPNAME = 'saucats'

# these variables are mandatory ('/' are converted automatically)
top, out = '.', 'build'



def options(context):
	context.load('compiler_cxx')



def configure(context):
	context.load('compiler_cxx')
	context.env.CXXFLAGS = ['-std=c++14', '-Wall', '-Wextra', '-O3', '-g', '-frounding-math']
	context.check_cfg(package = 'eigen3', uselib_store = 'eigen', args = ['--cflags'])
	context.check_cxx(lib = 'm', cflags = '-Wall', uselib_store = 'm')



def build(context):
	# saucats library
	context(
		name            = 'saucats',
		includes        = 'include',
		export_includes = 'include'
	)

	# unit-testing program
	context.program(
		target   = 'saucats-unitest',
		includes = 'tests',
		source   = context.path.ant_glob('tests/*.cpp'),
		use      = ['saucats', 'eigen', 'm']
	)

