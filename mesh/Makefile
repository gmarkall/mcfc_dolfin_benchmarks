all:
	./gmsh2triangle --2d cdisk.msh
	./dolfin-convert cdisk.msh cdisk.xml

remesh:
	gmsh -2 -clscale 0.32 -o cdisk.msh cdisk.geo

clean:
	rm -f *.ele *.edge *.node *.xml *.bak *.pyc

reallyclean: clean
	rf -f *.msh
