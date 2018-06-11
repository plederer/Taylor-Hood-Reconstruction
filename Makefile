NGSCXX = ngscxx
NGSLD = ngsld

objects = reconstruction_utils.o reconstruction_element.o reconstruction_vertex.o reconstruction_vertex_3D.o

headers = reconstruction.hpp

%.o : %.cpp $(headers)
	$(NGSCXX) -fdiagnostics-color=auto -Wno-unused-variable -I. -c $(subst .o,.cpp,$@) -o $@

reclib.so : $(objects)
	$(NGSLD) -shared $(objects) -lngsolve -lngcomp -lngfem -o $@

clean:
	rm *.o */*.o *.so 
