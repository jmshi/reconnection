#!/usr/bin/env python

from fipy import *
nx = 40
ny = nx
dx = 1./nx
dy = dx
L = dx * nx
mesh = Grid2D(dx=dx, dy=dy, nx=nx, ny=ny)

phi = CellVariable(name = "solution variable",
                   mesh = mesh,
                   value = 0.)
epsilon = 0.01
b = ((1,), (1,))
# (1) looks the same
# (2) looks a bit worse than (1)
# (3) different
convection = CentralDifferenceConvectionTerm # fails
convection = ExponentialConvectionTerm # (1)
convection = HybridConvectionTerm # (1)
convection = PowerLawConvectionTerm # (1)
#convection = UpwindConvectionTerm # (2)
#convection = ExplicitUpwindConvectionTerm # (3)
#convection = VanLeerConvectionTerm # (3)
eq = -convection(coeff=b)+DiffusionTerm(coeff=epsilon) == 0

x, y = mesh.getFaceCenters()
facesTopRight = ((mesh.getFacesRight() & (y > L / 2))
                | (mesh.getFacesTop() & (x > L / 2)))
facesBottomLeft = ((mesh.getFacesLeft() & (y < L / 2))
                    | (mesh.getFacesBottom() & (x < L / 2)))
BCs = (FixedValue(faces=facesTopRight, value=0),
       FixedValue(faces=facesBottomLeft, value=1))
eq.solve(var=phi, boundaryConditions=BCs)
print "1"
viewer = Viewer(vars=phi, datamin=0., datamax=1.)
#viewer.plot()
viewer.plotMesh()
viewer.plot(filename="plot1.png")
viewer.plot()
raw_input("Press <return> to proceed...")
