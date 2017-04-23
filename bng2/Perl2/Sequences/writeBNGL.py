from sequences import *
import os

a = NucObject.fromfile('seq.txt')
p = BNGNucObject.fromNucObject(a).__str__()


f = open('model.bngl','w')
strlist = []

strlist.append(
"""begin parameters
	k	1
	C	1
end parameters

begin molecule types
	d(t~A~T~C~G~X,p5,p3,c,b,fp~0~1)
	r(t~A~U~C~G~X,p5,p3,c,b,fp~0~1)
end molecule types

begin species"""
)

strlist.append("\t"+p+"\t1")
strlist.append(
"""end species

begin reaction rules"""
)

strlist.append("R1: "+p+"-> 0 k")
strlist.append(
"""end reaction rules

visualize({type=>"conventional"})"""
)


f.write("\n".join(strlist))
f.close()

os.system("perl C:\BioNetGenCodeBase\\bionetgen\\bng2\\BNG2.pl model.bngl")
