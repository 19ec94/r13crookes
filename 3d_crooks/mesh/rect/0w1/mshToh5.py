
#!/usr/bin/env python3

# pylint: disable=invalid-name

"""
Converter from geo-format to a mesh in h5-format.
Usage:
    python3 mesh_conversion_GeoToH5.py ring.geo ring0.h5 "-setnumber p 0"
    python3 mesh_conversion_GeoToH5.py ring.geo ring1.h5 "-setnumber p 1"
"""
import os
import sys
import dolfin as df

# Constants
GMSH_PATH = "gmsh"

if (len(sys.argv) < 3):
    print("Provide names of input.msh output.h5 files")
    exit(1)

# Read input/output file names
input_file_msh = sys.argv[1]
output_file_h5 = sys.argv[2]
intermediate_file_xml = input_file_msh.replace(".msh",".xml")  

# Create msh-mesh with Gmsh
#os.system( "{} {} -2 -o {}.msh {}".format(GMSH_PATH, gmsh_arguments, tmp_name, geo_input_file))

# Convert msh-format to xml-format using dolfin-convert program
os.system("dolfin-convert {0} {1}".format(input_file_msh,
    intermediate_file_xml )
    )

# Delete msh-mesh
#os.remove("{}.msh".format(tmp_name))

# Read xml-mesh
intermediate_file = intermediate_file_xml.replace(".xml", "")
mesh = df.Mesh("{}.xml".format(intermediate_file))
subdomains = df.MeshFunction(
        "size_t", mesh, 
        "{}_physical_region.xml".format( intermediate_file)
        )
boundaries = df.MeshFunction(
        "size_t", mesh, 
        "{}_facet_region.xml".format(intermediate_file)
        )
print("Created '.xml' files successfully")

# Delete xml-mesh
os.remove("{}.xml".format(intermediate_file))
os.remove("{}_physical_region.xml".format(intermediate_file))
os.remove("{}_facet_region.xml".format(intermediate_file))
print("Deleted '.xml' files successfully")

# Write h5-mesh
file = df.HDF5File(mesh.mpi_comm(), output_file_h5, "w")
file.write(mesh, "/mesh")
file.write(subdomains, "/subdomains")
file.write(boundaries, "/boundaries")
file.close()
