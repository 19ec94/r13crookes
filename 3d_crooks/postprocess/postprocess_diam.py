import sys
import  time as time_module
import yaml
import dolfin as df
import ufl
import tensoroperations as to
import mpi4py 

df.parameters['ghost_mode']= 'shared_facet' #for parallel program


# check if input file is given
if (len(sys.argv) != 2 ):
    print("provide an input yaml file")
    exit(1)

# read the input file
# TODO: check if the input file is a valid yaml file
yaml_file = sys.argv[1] 
with open (yaml_file, "r") as stream:
    inputfile= yaml.safe_load(stream)

# Check if the simulation results are present
output_folder = inputfile["output_folder"]+"/"
theta_in = df.XDMFFile("{0}theta.xdmf".format(output_folder))
s_in = df.XDMFFile("{0}s.xdmf".format(output_folder))
p_in = df.XDMFFile("{0}p.xdmf".format(output_folder))
u_in = df.XDMFFile("{0}u.xdmf".format(output_folder))
sigma_in = df.XDMFFile("{0}sigma.xdmf".format(output_folder))

# read the h5 mesh file name
h5_file = inputfile["meshes"][0]

# mesh - construct a dolfin mesh
sys.stdout.flush()
start_t = time_module.time()
mesh = df.Mesh()
hdf = df.HDF5File(mesh.mpi_comm(), h5_file, "r")
# mesh - read topology/subdomains/boundaries from mesh file
hdf.read(mesh, "/mesh", False)
dim = mesh.topology().dim()
subdomains = df.MeshFunction("size_t", mesh, dim)
hdf.read(subdomains, "/subdomains")
boundaries = df.MeshFunction("size_t", mesh, dim-1)
hdf.read(boundaries, "/boundaries")
end_t = time_module.time()
print("Finished reading mesh {}".format(end_t-start_t))
sys.stdout.flush()
print("*"*30)
print("Vertices ", mesh.num_vertices())
print("Cells ", mesh.num_cells())
print("hmin ", mesh.hmin())
print("hmax ", mesh.hmax())
print("*"*30)


df.dx = df.Measure("dx", domain=mesh, subdomain_data=subdomains)
df.ds = df.Measure("ds", domain=mesh, subdomain_data=boundaries)
df.dS = df.Measure("dS", domain=mesh, subdomain_data=boundaries)

# functionSpace - set up function spaces
cell = mesh.ufl_cell()
elems = {
        "theta": None,
        "s" : None,
        "p" : None,
        "u":None,
        "sigma" : None,
        "alpha" : None
        }
var_ranks = {
        "theta" : 0,
        "s" : 1,
        "p" : 0,
        "u" : 1,
        "sigma" : 2,
        "alpha" : 0
        }
mxd_elems = {
        "r13" : None
        }
fspaces = { 
           "theta" : None,
           "s" : None,
           "p" : None,
           "u" : None,
           "sigma" : None,
           "alpha" : None
           }
mxd_fspaces = {
        "r13" : None
        }

for var in elems:
    e = inputfile["elements"][var]["shape"]
    deg = inputfile["elements"][var]["degree"]
    if (var_ranks[var] == 0):
        elems[var] = df.FiniteElement(e, cell, deg)
    elif (var_ranks[var] == 1):
        elems[var] = df.VectorElement(e, cell, deg)
    elif (var_ranks[var] == 2):
        elems[var] = df.TensorElement(e, cell, deg,
#                                      symmetry ={
#                                          (0,1) : (1,0),
#                                          (2,0) : (0,2),
#                                          (1,2) : (2,1),
#                                          (2,2) : (0,0)
#                                          }
                                      )
    fspaces[var] = df.FunctionSpace(mesh, elems[var])

sys.stdout.flush()
start_t = time_module.time()
mxd_elems["r13"] = df.MixedElement(
        [ elems["theta"], elems["s"], 
         elems["p"], elems["u"], elems["sigma"],
         elems["alpha"]
         ]
        )
mxd_fspaces["r13"] = df.FunctionSpace(
        mesh,
        mxd_elems["r13"]
        )
end_t = time_module.time()
print("Finished creating function space {} ".format(end_t-start_t))
sys.stdout.flush()


sol  = df.Function(mxd_fspaces["r13"])
(theta,s,p,u,sigma,alpha) = sol.split(True)

# files are read above in variables. Copy that into normal fenics
# variables
theta_in.read_checkpoint(theta,"theta", 0)
s_in.read_checkpoint(s,"s", 0)
p_in.read_checkpoint(p,"p", 0)
u_in.read_checkpoint(u,"u", 0)
sigma_in.read_checkpoint(sigma,"sigma", 0)

#NOTE: The following projection is already done in FEniCSR13.
# symmetrise the stress tensor
#sys.stdout.flush()
#start_t = time_module.time()
#sigma = df.project(
#        to.gen3DTFdim3(sigma), df.TensorFunctionSpace(
#            mesh, "Lagrange",
#            inputfile["elements"]["sigma"]["degree"]
#            ), solver_type=inputfile["solver"]["solver_name"],
#        preconditioner_type=inputfile["solver"]["preconditioner"]
#        )
#end_t = time_module.time()
#sigma = to.gen3DTFdim3(sigma)
#print("Finished projecting the tensor {} ".format(end_t-start_t))
#sys.stdout.flush()

n_vec = df.FacetNormal(mesh)
# Get IDs of hot and cold sides
rh = 188 #right hot
rc = 183 #right cold
uh = 218 #up hot
uc = 213 #up cold
lh = 128 #left hot
lc = 123 #left cold
dh = 48 #down hot
dc = 43 #down cold
rtt = 182020 #right transition top
rtb = 18202 #right transition bottom
ltt = 122020 #left transition top
ltb = 12202 #left transition bottom
utt = 212020 #up transition top
utb = 21202 #up transition bottom
dtt = 42020 #down transition top
dtb = 4202 #down transition bottom
print("*"*30)
vane_surface_area =  df.assemble(1*df.ds(rh))
print("Integral of surface area is ", vane_surface_area )
print("Integral of pressure is", df.assemble(p*df.dx))

# Format flot
def ff(number):
    return "{:+0.4e}".format(number)

pressure_tensor = df.as_matrix([
    [sigma[0,0]+p, sigma[0,1],sigma[0,2]],
    [sigma[1,0], sigma[1,1]+p,sigma[1,2]],
    [sigma[2,0], sigma[2,1],sigma[2,2]+p],
    ])
pressure_tensor_proj = df.dot(pressure_tensor, n_vec)
def calculate_force_vector(surf_id):
    f_xx = df.assemble(pressure_tensor_proj[0]*df.ds(surf_id))
    f_yy = df.assemble(pressure_tensor_proj[1]*df.ds(surf_id))
    f_zz = df.assemble(pressure_tensor_proj[2]*df.ds(surf_id))
    return [f_xx, f_yy, f_zz]

# right vane
rh_f = calculate_force_vector(rh);
rc_f = calculate_force_vector(rc);
rtt_f = calculate_force_vector(rtt);
rtb_f = calculate_force_vector(rtb);
# up vane
uh_f = calculate_force_vector(uh);
uc_f = calculate_force_vector(uc);
utt_f = calculate_force_vector(utt);
utb_f = calculate_force_vector(utb);
# left vane
lh_f = calculate_force_vector(lh);
lc_f = calculate_force_vector(lc);
ltt_f = calculate_force_vector(ltt);
ltb_f = calculate_force_vector(ltb);
# down vane
dh_f = calculate_force_vector(dh);
dc_f = calculate_force_vector(dc);
dtt_f = calculate_force_vector(dtt);
dtb_f = calculate_force_vector(dtb);
#
print("*"*30)
print("Right vane")
print("RH : f_x, f_y, f_z: {}, {}, {}, Right Hot"\
        .format(ff(rh_f[0]), ff(rh_f[1]), ff(rh_f[2])))
print("RC : f_x, f_y, f_z: {}, {}, {}, Right Cold"\
        .format(ff(rc_f[0]), ff(rc_f[1]), ff(rc_f[2])))
print("NrF: f_x, f_y, f_z: {}, {}, {}, Normal Force"\
        .format(ff(rh_f[0]+rc_f[0]),
                ff(rh_f[1]+rc_f[1]),
                ff(rh_f[2]+rc_f[2])))
print("RTT: f_x, f_y, f_z: {}, {}, {}, Right Transition Top"\
        .format(ff(rtt_f[0]), ff(rtt_f[1]), ff(rtt_f[2])))
print("RTB: f_x, f_y, f_z: {}, {}, {}, Right Transition Bottom"\
        .format(ff(rtb_f[0]), ff(rtb_f[1]), ff(rtb_f[2])))
print("TrF: f_x, f_y, f_z: {}, {}, {}, Transition Force"\
        .format(ff(rtt_f[0]+rtb_f[0]), 
                ff(rtt_f[1]+rtb_f[1]),
                ff(rtt_f[2]+rtb_f[2])))
print("NF : f_x, f_y, f_z: {}, {}, {}, Net Force"\
        .format(ff(rh_f[0]+rc_f[0]+rtt_f[0]+rtb_f[0]),
                ff(rh_f[1]+rc_f[1]+rtt_f[1]+rtb_f[1]),
                ff(rh_f[2]+rc_f[2]+rtt_f[2]+rtb_f[2])))
print("Left vane")
print("RH : f_x, f_y, f_z: {}, {}, {}, Left Hot"\
        .format(ff(lh_f[0]), ff(lh_f[1]), ff(lh_f[2])))
print("LC : f_x, f_y, f_z: {}, {}, {}, Left Cold"\
        .format(ff(lc_f[0]), ff(lc_f[1]), ff(lc_f[2])))
print("NrF: f_x, f_y, f_z: {}, {}, {}, Normal Force"\
        .format(ff(lh_f[0]+lc_f[0]), 
                ff(lh_f[1]+lc_f[1]),
                ff(lh_f[2]+lc_f[2])))
print("LTT: f_x, f_y, f_z: {}, {}, {}, Left Transition Top"\
        .format(ff(ltt_f[0]), ff(ltt_f[1]), ff(ltt_f[2])))
print("LTB: f_x, f_y, f_z: {}, {}, {}, Left Transition Bottom"\
        .format(ff(ltb_f[0]), ff(ltb_f[1]), ff(ltb_f[2])))
print("TrF: f_x, f_y, f_z: {}, {}, {}, Transition Force"\
        .format(ff(ltt_f[0]+ltb_f[0]), 
                ff(ltt_f[1]+ltb_f[1]),
                ff(ltt_f[2]+ltb_f[2])))
print("NF : f_x, f_y, f_z: {}, {}, {}, Net Force"\
        .format(ff(lh_f[0]+lc_f[0]+ltt_f[0]+ltb_f[0]), 
                ff(lh_f[1]+lc_f[1]+ltt_f[1]+ltb_f[1]),
                ff(lh_f[2]+lc_f[2]+ltt_f[2]+ltb_f[2])))
print("Up vane")
print("UH : f_x, f_y, f_z: {}, {}, {}, Up Hot"\
        .format(ff(uh_f[0]), ff(uh_f[1]), ff(uh_f[2])))
print("UC : f_x, f_y, f_z: {}, {}, {}, Up Cold"\
        .format(ff(uc_f[0]), ff(uc_f[1]), ff(uc_f[2])))
print("NrF: f_x, f_y, f_z: {}, {}, {}, Normal Force"\
        .format(ff(uh_f[0]+uc_f[0]),
                ff(uh_f[1]+uc_f[1]),
                ff(uh_f[2]+uc_f[2])))
print("UTT: f_x, f_y, f_z: {}, {}, {}, Up Transition Top"\
        .format(ff(utt_f[0]), ff(utt_f[1]), ff(utt_f[2])))
print("UTB: f_x, f_y, f_z: {}, {}, {}, Up Transition Bottom"\
        .format(ff(utb_f[0]), ff(utb_f[1]), ff(utb_f[2])))
print("TrF: f_x, f_y, f_z: {}, {}, {}, Transition Force"\
        .format(ff(utt_f[0]+utb_f[0]), 
                ff(utt_f[1]+utb_f[1]),
                ff(utt_f[2]+utb_f[2])))
print("NF : f_x, f_y, f_z: {}, {}, {}, Net Force"\
        .format(ff(uh_f[0]+uc_f[0]+utt_f[0]+utb_f[0]),
                ff(uh_f[1]+uc_f[1]+utt_f[1]+utb_f[1]),
                ff(uh_f[2]+uc_f[2]+utt_f[2]+utb_f[2])))
print("Down vane")
print("DH : f_x, f_y, f_z: {}, {}, {}, Down Hot"\
        .format(ff(dh_f[0]), ff(dh_f[1]), ff(dh_f[2])))
print("DC : f_x, f_y, f_z: {}, {}, {}, Down Cold"\
        .format(ff(dc_f[0]), ff(dc_f[1]), ff(dc_f[2])))
print("NrF: f_x, f_y, f_z: {}, {}, {}, Normal Force"\
        .format(ff(dh_f[0]+dc_f[0]),
                ff(dh_f[1]+dc_f[1]),
                ff(dh_f[2]+dc_f[2])))
print("DTT: f_x, f_y, f_z: {}, {}, {}, Down Transition Top"\
        .format(ff(dtt_f[0]), ff(dtt_f[1]), ff(dtt_f[2])))
print("DTB: f_x, f_y, f_z: {}, {}, {}, Down Transition Bottom"\
        .format(ff(dtb_f[0]), ff(dtb_f[1]), ff(dtb_f[2])))
print("TrF: f_x, f_y, f_z: {}, {}, {}, Transition Force"\
        .format(ff(dtt_f[0]+dtb_f[0]), 
                ff(dtt_f[1]+dtb_f[1]),
                ff(dtt_f[2]+dtb_f[2])))
print("NF : f_x, f_y, f_z: {}, {}, {}, Net Force"\
        .format(ff(dh_f[0]+dc_f[0]+dtt_f[0]+dtb_f[0]),
                ff(dh_f[1]+dc_f[1]+dtt_f[1]+dtb_f[1]),
                ff(dh_f[2]+dc_f[2]+dtt_f[2]+dtb_f[2])))
print("*"*30)
with open("{0}force_recalculation.txt".format(output_folder), "w") as forcefile:
    forcefile.write("Right vane\n")
    forcefile.write("RH : f_x, f_y, f_z: {}, {}, {}, Right Hot\n"\
            .format(ff(rh_f[0]), ff(rh_f[1]), ff(rh_f[2])))
    forcefile.write("RC : f_x, f_y, f_z: {}, {}, {}, Right Cold\n"\
            .format(ff(rc_f[0]), ff(rc_f[1]), ff(rc_f[2])))
    forcefile.write("NrF: f_x, f_y, f_z: {}, {}, {}, Normal Force\n"\
            .format(ff(rh_f[0]+rc_f[0]),
                    ff(rh_f[1]+rc_f[1]),
                    ff(rh_f[2]+rc_f[2])))
    forcefile.write("RTT: f_x, f_y, f_z: {}, {}, {}, Right Transition Top\n"\
            .format(ff(rtt_f[0]), ff(rtt_f[1]), ff(rtt_f[2])))
    forcefile.write("RTB: f_x, f_y, f_z: {}, {}, {}, Right Transition Bottom\n"\
            .format(ff(rtb_f[0]), ff(rtb_f[1]), ff(rtb_f[2])))
    forcefile.write("TrF: f_x, f_y, f_z: {}, {}, {}, Transition Force\n"\
            .format(ff(rtt_f[0]+rtb_f[0]), 
                    ff(rtt_f[1]+rtb_f[1]),
                    ff(rtt_f[2]+rtb_f[2])))
    forcefile.write("NF : f_x, f_y, f_z: {}, {}, {}, Net Force\n"\
            .format(ff(rh_f[0]+rc_f[0]+rtt_f[0]+rtb_f[0]),
                    ff(rh_f[1]+rc_f[1]+rtt_f[1]+rtb_f[1]),
                    ff(rh_f[2]+rc_f[2]+rtt_f[2]+rtb_f[2])))
    forcefile.write("Left vane\n")
    forcefile.write("LH : f_x, f_y, f_z: {}, {}, {}, Left Hot\n"\
            .format(ff(lh_f[0]), ff(lh_f[1]), ff(lh_f[2])))
    forcefile.write("LC : f_x, f_y, f_z: {}, {}, {}, Left Cold\n"\
            .format(ff(lc_f[0]), ff(lc_f[1]), ff(lc_f[2])))
    forcefile.write("NrF: f_x, f_y, f_z: {}, {}, {}, Normal Force\n"\
            .format(ff(lh_f[0]+lc_f[0]), 
                    ff(lh_f[1]+lc_f[1]),
                    ff(lh_f[2]+lc_f[2])))
    forcefile.write("LTT: f_x, f_y, f_z: {}, {}, {}, Left Transition Top\n"\
            .format(ff(ltt_f[0]), ff(ltt_f[1]), ff(ltt_f[2])))
    forcefile.write("LTB: f_x, f_y, f_z: {}, {}, {}, Left Transition Bottom\n"\
            .format(ff(ltb_f[0]), ff(ltb_f[1]), ff(ltb_f[2])))
    forcefile.write("TrF: f_x, f_y, f_z: {}, {}, {}, Transition Force\n"\
            .format(ff(ltt_f[0]+ltb_f[0]), 
                    ff(ltt_f[1]+ltb_f[1]),
                    ff(ltt_f[2]+ltb_f[2])))
    forcefile.write("NF : f_x, f_y, f_z: {}, {}, {}, Net Force\n"\
            .format(ff(lh_f[0]+lc_f[0]+ltt_f[0]+ltb_f[0]), 
                    ff(lh_f[1]+lc_f[1]+ltt_f[1]+ltb_f[1]),
                    ff(lh_f[2]+lc_f[2]+ltt_f[2]+ltb_f[2])))
    forcefile.write("Up vane\n")
    forcefile.write("UH : f_x, f_y, f_z: {}, {}, {}, Up Hot\n"\
            .format(ff(uh_f[0]), ff(uh_f[1]), ff(uh_f[2])))
    forcefile.write("UC : f_x, f_y, f_z: {}, {}, {}, Up Cold\n"\
            .format(ff(uc_f[0]), ff(uc_f[1]), ff(uc_f[2])))
    forcefile.write("NrF: f_x, f_y, f_z: {}, {}, {}, Normal Force\n"\
            .format(ff(uh_f[0]+uc_f[0]),
                    ff(uh_f[1]+uc_f[1]),
                    ff(uh_f[2]+uc_f[2])))
    forcefile.write("UTT: f_x, f_y, f_z: {}, {}, {}, Up Transition Top\n"\
            .format(ff(utt_f[0]), ff(utt_f[1]), ff(utt_f[2])))
    forcefile.write("UTB: f_x, f_y, f_z: {}, {}, {}, Up Transition Bottom\n"\
            .format(ff(utb_f[0]), ff(utb_f[1]), ff(utb_f[2])))
    forcefile.write("TrF: f_x, f_y, f_z: {}, {}, {}, Transition Force\n"\
            .format(ff(utt_f[0]+utb_f[0]),
                    ff(utt_f[1]+utb_f[1]),
                    ff(utt_f[2]+utb_f[2])))
    forcefile.write("NF : f_x, f_y, f_z: {}, {}, {}, Net Force\n"\
            .format(ff(uh_f[0]+uc_f[0]+utt_f[0]+utb_f[0]),
                    ff(uh_f[1]+uc_f[1]+utt_f[1]+utb_f[1]),
                    ff(uh_f[2]+uc_f[2]+utt_f[2]+utb_f[2])))
    forcefile.write("Down vane\n")
    forcefile.write("DH : f_x, f_y, f_z: {}, {}, {}, Down Hot side\n"\
            .format(ff(dh_f[0]), ff(dh_f[1]), ff(dh_f[2])))
    forcefile.write("DC : f_x, f_y, f_z: {}, {}, {}, Down Cold side\n"\
            .format(ff(dc_f[0]), ff(dc_f[1]), ff(dc_f[2])))
    forcefile.write("NrF: f_x, f_y, f_z: {}, {}, {}, Normal Force\n"\
            .format(ff(dh_f[0]+dc_f[0]),
                    ff(dh_f[1]+dc_f[1]),
                    ff(dh_f[2]+dc_f[2])))
    forcefile.write("DTT: f_x, f_y, f_z: {}, {}, {}, Down Transition Top\n"\
            .format(ff(dtt_f[0]), ff(dtt_f[1]), ff(dtt_f[2])))
    forcefile.write("DTB: f_x, f_y, f_z: {}, {}, {}, Down Transition Bottom\n"\
            .format(ff(dtb_f[0]), ff(dtb_f[1]), ff(dtb_f[2])))
    forcefile.write("TrF: f_x, f_y, f_z: {}, {}, {}, Transition Force\n"\
            .format(ff(dtt_f[0]+dtb_f[0]),
                    ff(dtt_f[1]+dtb_f[1]),
                    ff(dtt_f[2]+dtb_f[2])))
    forcefile.write("NF : f_x, f_y, f_z: {}, {}, {}, Net Force\n"\
            .format(ff(dh_f[0]+dc_f[0]+dtt_f[0]+dtb_f[0]),
                    ff(dh_f[1]+dc_f[1]+dtt_f[1]+dtb_f[1]),
                    ff(dh_f[2]+dc_f[2]+dtt_f[2]+dtb_f[2])))

print("Successfully written forces to file force.txt")
