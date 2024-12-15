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
        }
var_ranks = {
        "theta" : 0,
        "s" : 1,
        "p" : 0,
        "u" : 1,
        "sigma" : 2,
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
                                      )
    fspaces[var] = df.FunctionSpace(mesh, elems[var])

sys.stdout.flush()
start_t = time_module.time()
mxd_elems["r13"] = df.MixedElement(
        [ elems["theta"], elems["s"], 
         elems["p"], elems["u"], elems["sigma"]
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
(theta,s,p,u,sigma) = sol.split(True)

# files are read above in variables. Copy that into normal fenics
# variables
#theta_in.read_checkpoint(theta,"theta", 0)
#s_in.read_checkpoint(s,"s", 0)
p_in.read_checkpoint(p,"p", 0)
#u_in.read_checkpoint(u,"u", 0)
sigma_in.read_checkpoint(sigma,"sigma", 0)

print("*"*30)
n_vec = df.FacetNormal(mesh)
# Get IDs of hot and cold sides
#
rh = 188 #right hot
uh = 218 #up hot
lh = 128 #left hot
dh = 48 #down hot
#
rc = 183 #right cold
uc = 213 #up cold
lc = 123 #left cold
dc = 43 #down cold
#
rt = 1820 #right transition 
lt = 1220 #left transition
ut = 2120 # up transition
dt = 420 #down trnasition

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
rt_f = calculate_force_vector(rt);
# up vane
uh_f = calculate_force_vector(uh);
uc_f = calculate_force_vector(uc);
ut_f = calculate_force_vector(ut);
# left vane
lh_f = calculate_force_vector(lh);
lc_f = calculate_force_vector(lc);
lt_f = calculate_force_vector(lt);
# down vane
dh_f = calculate_force_vector(dh);
dc_f = calculate_force_vector(dc);
dt_f = calculate_force_vector(dt);

print("Right vane")
print("RH : f_x, f_y, f_z: {}, {}, {}, Right Hot"\
        .format(ff(rh_f[0]), ff(rh_f[1]), ff(rh_f[2])))
print("RC : f_x, f_y, f_z: {}, {}, {}, Right Cold"\
        .format(ff(rc_f[0]), ff(rc_f[1]), ff(rc_f[2])))
print("NrF: f_x, f_y, f_z: {}, {}, {}, Normal Force"\
        .format(ff(rh_f[0]+rc_f[0]),
                ff(rh_f[1]+rc_f[1]),
                ff(rh_f[2]+rc_f[2])))
print("RT : f_x, f_y, f_z: {}, {}, {}, Right Transition"\
        .format(ff(rt_f[0]), ff(rt_f[1]), ff(rt_f[2])))
print("NF : f_x, f_y, f_z: {}, {}, {}, Net Force"\
        .format(ff(rh_f[0]+rc_f[0]+rt_f[0]),
                ff(rh_f[1]+rc_f[1]+rt_f[1]),
                ff(rh_f[2]+rc_f[2]+rt_f[2])))
print("Left vane")
print("RH : f_x, f_y, f_z: {}, {}, {}, Left Hot"\
        .format(ff(lh_f[0]), ff(lh_f[1]), ff(lh_f[2])))
print("LC : f_x, f_y, f_z: {}, {}, {}, Left Cold"\
        .format(ff(lc_f[0]), ff(lc_f[1]), ff(lc_f[2])))
print("NrF : f_x, f_y, f_z: {}, {}, {}, Normal Force"\
        .format(ff(lh_f[0]+lc_f[0]),
                ff(lh_f[1]+lc_f[1]),
                ff(lh_f[2]+lc_f[2])))
print("LT : f_x, f_y, f_z: {}, {}, {}, Left Transition"\
        .format(ff(lt_f[0]), ff(lt_f[1]), ff(lt_f[2])))
print("NF : f_x, f_y, f_z: {}, {}, {}, Net Force"\
        .format(ff(lh_f[0]+lc_f[0]+lt_f[0]),
                ff(lh_f[1]+lc_f[1]+lt_f[1]),
                ff(lh_f[2]+lc_f[2]+lt_f[2])))
print("Up vane")
print("UH : f_x, f_y, f_z: {}, {}, {}, Up Hot"\
        .format(ff(uh_f[0]), ff(uh_f[1]), ff(uh_f[2])))
print("UC : f_x, f_y, f_z: {}, {}, {}, Up Cold"\
        .format(ff(uc_f[0]), ff(uc_f[1]), ff(uc_f[2])))
print("NrF: f_x, f_y, f_z: {}, {}, {}, Normal Force"\
        .format(ff(uh_f[0]+uc_f[0]),
                ff(uh_f[1]+uc_f[1]),
                ff(uh_f[2]+uc_f[2])))
print("UT : f_x, f_y, f_z: {}, {}, {}, Up Transition"\
        .format(ff(ut_f[0]), ff(ut_f[1]), ff(ut_f[2])))
print("NF : f_x, f_y, f_z: {}, {}, {}, Net Force"\
        .format(ff(uh_f[0]+uc_f[0]+ut_f[0]),
                ff(uh_f[1]+uc_f[1]+ut_f[1]),
                ff(uh_f[2]+uc_f[2]+ut_f[2])))
print("Down vane")
print("DH : f_x, f_y, f_z: {}, {}, {}, Down Hot"\
        .format(ff(dh_f[0]), ff(dh_f[1]), ff(dh_f[2])))
print("DC : f_x, f_y, f_z: {}, {}, {}, Down Cold"\
        .format(ff(dc_f[0]), ff(dc_f[1]), ff(dc_f[2])))
print("NFr: f_x, f_y, f_z: {}, {}, {}, Normal Force"\
        .format(ff(dh_f[0]+dc_f[0]),
                ff(dh_f[1]+dc_f[1]),
                ff(dh_f[2]+dc_f[2])))
print("DT : f_x, f_y, f_z: {}, {}, {}, Down Transition"\
        .format(ff(dt_f[0]), ff(dt_f[1]), ff(dt_f[2])))
print("NF : f_x, f_y, f_z: {}, {}, {}, Net Force"\
        .format(ff(dh_f[0]+dc_f[0]+dt_f[0]),
                ff(dh_f[1]+dc_f[1]+dt_f[1]),
                ff(dh_f[2]+dc_f[2]+dt_f[2])))
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
    forcefile.write("RT : f_x, f_y, f_z: {}, {}, {}, Right Transition\n"\
            .format(ff(rt_f[0]), ff(rt_f[1]), ff(rt_f[2])))
    forcefile.write("NF : f_x, f_y, f_z: {}, {}, {}, Net Force\n"\
            .format(ff(rh_f[0]+rc_f[0]+rt_f[0]),
                    ff(rh_f[1]+rc_f[1]+rt_f[1]),
                    ff(rh_f[2]+rc_f[2]+rt_f[2])))
    forcefile.write("Left vane\n")
    forcefile.write("LH : f_x, f_y, f_z: {}, {}, {}, Left Hot\n"\
            .format(ff(lh_f[0]), ff(lh_f[1]), ff(lh_f[2])))
    forcefile.write("LC : f_x, f_y, f_z: {}, {}, {}, Left Cold\n"\
            .format(ff(lc_f[0]), ff(lc_f[1]), ff(lc_f[2])))
    forcefile.write("NFr: f_x, f_y, f_z: {}, {}, {}, Normal Force\n"\
            .format(ff(lh_f[0]+lc_f[0]),
                    ff(lh_f[1]+lc_f[1]),
                    ff(lh_f[2]+lc_f[2])))
    forcefile.write("LT : f_x, f_y, f_z: {}, {}, {}, Left Transition\n"\
            .format(ff(lt_f[0]), ff(lt_f[1]), ff(lt_f[2])))
    forcefile.write("NF : f_x, f_y, f_z: {}, {}, {}, Net Force\n"\
            .format(ff(lh_f[0]+lc_f[0]+lt_f[0]),
                    ff(lh_f[1]+lc_f[1]+lt_f[1]),
                    ff(lh_f[2]+lc_f[2]+lt_f[2])))
    forcefile.write("Up vane\n")
    forcefile.write("UH : f_x, f_y, f_z: {}, {}, {}, Up Hot\n"\
            .format(ff(uh_f[0]), ff(uh_f[1]), ff(uh_f[2])))
    forcefile.write("UC : f_x, f_y, f_z: {}, {}, {}, Up Cold\n"\
            .format(ff(uc_f[0]), ff(uc_f[1]), ff(uc_f[2])))
    forcefile.write("NFr: f_x, f_y, f_z: {}, {}, {}, Normal Force\n"\
            .format(ff(uh_f[0]+uc_f[0]),
                    ff(uh_f[1]+uc_f[1]),
                    ff(uh_f[2]+uc_f[2])))
    forcefile.write("UT : f_x, f_y, f_z: {}, {}, {}, Up Transition\n"\
            .format(ff(ut_f[0]), ff(ut_f[1]), ff(ut_f[2])))
    forcefile.write("NF : f_x, f_y, f_z: {}, {}, {}, Net Force\n"\
            .format(ff(uh_f[0]+uc_f[0]+ut_f[0]),
                    ff(uh_f[1]+uc_f[1]+ut_f[1]),
                    ff(uh_f[2]+uc_f[2]+ut_f[2])))
    forcefile.write("Down vane\n")
    forcefile.write("DH : f_x, f_y, f_z: {}, {}, {}, Down Hot\n"\
            .format(ff(dh_f[0]), ff(dh_f[1]), ff(dh_f[2])))
    forcefile.write("DC : f_x, f_y, f_z: {}, {}, {}, Down Cold\n"\
            .format(ff(dc_f[0]), ff(dc_f[1]), ff(dc_f[2])))
    forcefile.write("NrF: f_x, f_y, f_z: {}, {}, {}, Normal Force\n"\
            .format(ff(dh_f[0]+dc_f[0]),
                    ff(dh_f[1]+dc_f[1]),
                    ff(dh_f[2]+dc_f[2])))
    forcefile.write("DT : f_x, f_y, f_z: {}, {}, {}, Down Transition\n"\
            .format(ff(dt_f[0]), ff(dt_f[1]), ff(dt_f[2])))
    forcefile.write("NF : f_x, f_y, f_z: {}, {}, {}, Net Force\n"\
            .format(ff(dh_f[0]+dc_f[0]+dt_f[0]),
                    ff(dh_f[1]+dc_f[1]+dt_f[1]),
                    ff(dh_f[2]+dc_f[2]+dt_f[2])))

