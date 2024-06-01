
CUT_Y = 1000
CUT_Z = 3

for iy in range(CUT_Y):
    for iz in range(CUT_Z):
        print ("./Roessler_DDE_proof_piece out=grid3 src=grid3 dst=grid3 cuty={} cutz={} iy={} iz={}".format(CUT_Y, CUT_Z, iy, iz))


CUT_Y = 10
CUT_Z = 3
for iy in range(CUT_Y):
    for iz in range(CUT_Z):
        print ("./Roessler_DDE_proof_piece out=Pc3_0 'src=c3[0]' 'dst=c3[1]' cuty={} cutz={} iy={} iz={}".format(CUT_Y, CUT_Z, iy, iz))

CUT_Y = 40
CUT_Z = 3
for iy in range(CUT_Y):
    for iz in range(CUT_Z):
        print ("./Roessler_DDE_proof_piece out=Pc3_1 'src=c3[1]' 'dst=c3[2]' cuty={} cutz={} iy={} iz={}".format(CUT_Y, CUT_Z, iy, iz))

CUT_Y = 20
CUT_Z = 3
for iy in range(CUT_Y):
    for iz in range(CUT_Z):
        print ("./Roessler_DDE_proof_piece out=Pc3_2 'src=c3[2]' 'dst=c3[0]' cuty={} cutz={} iy={} iz={}".format(CUT_Y, CUT_Z, iy, iz))