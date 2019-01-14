#! /usr/bin/env python

# the following two variables are used by the target "waf dist"
VERSION = '1.0.0'
APPNAME = 'saucats'

# these variables are mandatory ('/' are converted automatically)
top, out = '.', 'build'



def options(context):
	context.load('compiler_cxx compiler_c gnu_dirs')



def configure(context):
	context.setenv('debug')
	context.load('compiler_cxx')
	context.env.CXXFLAGS = ['-std=c++14', '-Wall', '-Wextra', '-O3', '-g', '-frounding-math']
	context.check_cfg(package = 'eigen3', uselib_store = 'eigen', args = ['--cflags'])
	context.check_cc(lib = 'm', cflags = '-Wall', uselib_store = 'm')

	context.setenv('release', env = context.env.derive())



def init(context):
	from waflib.Build import BuildContext, CleanContext, InstallContext, UninstallContext

	for x in 'debug release'.split():
		for y in (BuildContext, CleanContext, InstallContext, UninstallContext):
			name = y.__name__.replace('Context','').lower()
			class tmp(y):
				cmd = name + '_' + x
				variant = x



def build(context):
	if not context.variant:
		context.fatal('Call "waf build_debug" or "waf build_release", and read the comments in the wscript file!')

	# saucats library
	context(
		name            = 'saucats',
		includes        = 'include',
		export_includes = 'include'
	)

	context.install_files(
		'${PREFIX}',
		context.path.ant_glob('include/saucats/**/*.h') + ['include/saucats/Geometry', 'include/saucats/Macros', 'include/saucats/SDF'],
		relative_trick = True
	)

	# unit-testing program
	if context.variant == 'debug':
		context.program(
			target   = 'saucats-unitest',
			includes = 'tests',
			source   = context.path.ant_glob('tests/*.cpp'),
			use      = ['saucats', 'eigen', 'm']
		)

