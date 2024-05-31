# this will make it into a separate subfolder of a given name, trailing / is important!) 
PACKAGE =  sharkovski-roessler-dde/

# the defaut output folder (defined above) is ignored if those two are given (please keep trailing '/' !)
OUTDIR = ./bin/
OBJDIR = ./.obj/

# a list of all the programs to compile (from this directory), separated by single space.
# i.e. .cpp files that contain main() function.
EXECUTABLES = Roessler_DDE_gencoords Roessler_DDE_prepR Roessler_DDE_prepPiece Roessler_DDE_proof Roessler_DDE_proof_piece

# a list of all your units to be linked with your programs (other .cpp files), separated by single space
# i.e. all other .cpp files that does not have main() inside.
OTHERS = utils

# 'run' the true Makefile
include ../makefile-run.mk
