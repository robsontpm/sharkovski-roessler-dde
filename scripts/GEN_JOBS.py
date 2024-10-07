# this file is to generate list of jobs to be run to
# generate the tail of the section and later to prove
# the covering relations. We do it that way to be able
# to send single tasks to separate processors on the
# level of the operating system and to be able to redo 
# some computations later, without affecting results from 
# other. This is becouse DDE computations tend to be very long
# and we should have a way to restart from where we stopped,
# or redo just a fraction of them.

# all folders are with respect to the ./ folder (scripts)
# unless otherwise stated.
DATA_DIR="../proof-data/"       # here will be the initial data for the proofs stored
                                # the PREParation scripts will put final data here
                                # the PROOF scripts will read data from here
TMP_DIR="../tmp/"               # the tmp forlder for PREParation scripts
PROOF_DIR="../proof-output/"    # common working directory, all data of the proof will go here
SCRIPT_DIR="../bin/"            # script directory w.r.t. to WORKING_DIRECTORY
# this is a common setup for jobs to prove the covering.
# in short, each entry { ... } defines the covering src => dst
# that should be proven with given number of subdivisions in y and z
# directions on the section x=0 
# (this subdivision applies to the head of the DDE solution! See paper.)
# here we have 4 coverings: G3=>G3, Ci=>C{i+1} (modulo 3)
# the output of the programs should go to 'out' directory (this is technical)
COVERING_SETUP=[
    { "src": "G", "dst": "G", "out": "G", "cut_y": 1000, "cut_z": 3, },
    { "src": "S_0", "dst": "S_1", "out": "PS_1", "cut_y": 10, "cut_z": 3, },
    { "src": "S_1", "dst": "S_2", "out": "PS_1", "cut_y": 40, "cut_z": 3, },
    { "src": "S_2", "dst": "S_0", "out": "PS_2", "cut_y": 20, "cut_z": 3, },
]
# now we setup two major jobs: first is to prepare set,
# the second job is to just run the proof. 
JOBS = [
    {
        "filepath": "./PREP_JOBLIST.sh",
        "script": SCRIPT_DIR + "Roessler_DDE_prepare_piece",
        "wd": TMP_DIR,
        "outdir": DATA_DIR,
        "subjobs": COVERING_SETUP,
    },
    {
        "filepath": "./PROOF_JOBLIST.sh",
        "script": SCRIPT_DIR + "Roessler_DDE_proof_piece",
        "wd": DATA_DIR,
        "outdir": PROOF_DIR,        
        "subjobs": COVERING_SETUP,
    },
]
# please do not modify this
COMMAND="{script} 'wd={wd}' 'out={outdir}{out}' 'src={src}' 'dst={dst}' 'cuty={cut_y}' 'cutz={cut_z}' 'iy={iy}' 'iz={iz}'\n"

for setup in JOBS:    
    # update local variables filepath, script, subjobs
    # this is considered unsafe, but everybody will be doing it locally on their own computer
    locals().update(setup) 
    outf = open(filepath, "w")
    for item in subjobs:
        sjob = dict(item)       # make a full copy
        sjob.update(setup)      # insert script into job        
        for iy in range(sjob["cut_y"]):
            for iz in range(sjob["cut_z"]):
                sjob.update({ "iy": iy, "iz": iz, })
                outf.write(COMMAND.format(**sjob))
    outf.close()
