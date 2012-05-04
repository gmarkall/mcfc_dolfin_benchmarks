all:
	gmsh -2 -clscale 0.32 -o cdisk.msh cdisk.geo
	gmsh2triangle --2d cdisk.msh
	./dolfin-convert cdisk.msh cdisk.xml

clean:
	rm -f *.ele *.edge *.node *.msh *.xml *.bak *.pyc
