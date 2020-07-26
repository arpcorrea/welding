from consolidate import *



cwd = os.getcwd()
# All the data from deck.yaml is now in the following deck variable

deck = Deck( cwd + "/TwoPlates.yaml" )

geometry = Geometry(deck)

meshes = MeshTwoPlates( deck,geometry )

model_HT= HeatTransfer(deck,meshes)

model_IC = IntimateContact(meshes,deck)

plots=PlotsTwoPlates(deck,meshes,meshes.T,meshes.Dic)

solves = SolvesTwoPlates( deck,model_HT,meshes,plots,model_IC)


