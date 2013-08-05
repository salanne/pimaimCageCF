#env = Environment(F90='mpif90',LINK='mpif90',LINKFLAGS='-Ofast -funroll-loops -march=native',F90FLAGS='-Ofast -funroll-loops -march=native')  # Initialize the environment
env = Environment(F90='mpif90', LINK='mpif90', LINKFLAGS='', F90FLAGS='-Jobj -Wall')
# The next line of code is an array of the source files names used in the program.
sources = Glob('src/*.f90')

# The next line is the actual code that links the executable. env.Program generates an executable

objs = env.Program('pimaimCageCF', sources)
