from numpy.distutils.fcompiler import new_fcompiler

compiler = new_fcompiler(compiler='intelem')
compiler.customize()
compiler.show_customization()

