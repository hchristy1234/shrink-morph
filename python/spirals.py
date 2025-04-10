import shrink_morph_py
import numpy as np
import polyscope as ps
import polyscope.imgui as gui
import togcode
import tkinter
from tkinter import filedialog
import triangle as tr

root = tkinter.Tk()
root.withdraw()

def display_trajectories(nodes, edges):
  ps_traj = ps.register_curve_network("Layer", nodes, edges, enabled=True, radius=printer.nozzle_width / 2)
  ps_traj.set_radius(printer.nozzle_width / 2, relative=False)
  ps_traj.set_color((0.3, 0.6, 0.8))


def convert_trajectories(trajectories):
  nodes = trajectories[0]
  edges = np.empty([nodes.shape[0] - 1, 2])
  edges[:, 0] = np.arange(nodes.shape[0] - 1)
  edges[:, 1] = np.arange(1, nodes.shape[0])

  for traj in trajectories:
    new_edges = np.empty([traj.shape[0] - 1, 2])
    new_edges[:, 0] = np.arange(nodes.shape[0], nodes.shape[0] + traj.shape[0] - 1)
    new_edges[:, 1] = np.arange(nodes.shape[0] + 1, nodes.shape[0] + traj.shape[0])
    
    nodes = np.vstack((nodes, traj))
    edges = np.vstack((edges, new_edges))

  return nodes, edges

printer_profile = "Prusa_MK3S"
printer = togcode.Printer(printer_profile)
printer.layer_height = 0.2

lambda1 = 1
lambda2 = 0.92
outer_radius = 30
eccentricity = 1
inner_radius = 1
thickness = 0.4
angle = 45

vertices = []
segments = []
holes = [[0,0]]

def make_circle(radius, factor):
  n = len(vertices)
  n_circle = round(factor * radius * np.pi)
  for i in range(n_circle):
    vertices.append([radius * np.cos(2 * i / n_circle * np.pi), radius * np.sin(2 * i / n_circle * np.pi)])
    segments.append([n + i, n + ((i + 1) % n_circle)])

make_circle(outer_radius, 3)
make_circle(inner_radius, 3) 

A = dict(vertices=vertices, segments=segments, holes=holes)
B = tr.triangulate(A, 'pqa0.1')

# Initialize polyscope
ps.init()
ps.set_give_focus_on_show(True)
ps.set_ground_plane_mode("shadow_only")

zeros = np.zeros((len(B['vertices']), 1))
P = np.hstack((np.array(B['vertices']), zeros))
F = np.array(B['triangles'])
ps.register_surface_mesh("Parameterization", P, F, edge_width=1, color=(42/255, 53/255, 213/255))

theta = np.empty(P.shape[0])
def callback():
  global thickness, outer_radius, eccentricity, angle, lambda2, P, trajectories, printer, vertices, segments, F

  gui.PushItemWidth(50)
  _, printer.nozzle_width = gui.InputDouble("Nozzle width (mm)", printer.nozzle_width, format="%.2f")    
  _, printer.print_speed = gui.InputDouble("Printing speed (mm/s)", printer.print_speed, format="%.0f")    
  _, printer.first_layer_speed = gui.InputDouble("First layer speed (mm/s)", printer.first_layer_speed, format="%.0f")
  _, thickness = gui.InputDouble("Total thickness", thickness, format="%.2f")
  _, printer.layer_height = gui.InputDouble("Layer height", printer.layer_height, format="%.2f")
  _, printer.bed_temp = gui.InputDouble("Bed temperature (ºC)", printer.bed_temp, format="%.0f")
  _, printer.extruder_temp = gui.InputDouble("Nozzle temperature (ºC)", printer.extruder_temp, format="%.0f")
  _, printer.filament_priming = gui.InputDouble("Retract length (mm)", printer.filament_priming, format="%.1f")
  gui.PushItemWidth(110)
  _, printer.nloops = gui.InputInt("Number of loops around object", printer.nloops, step=1)
  gui.PopItemWidth()
  _, printer.flow_multiplier = gui.InputDouble("Flow multiplier", printer.flow_multiplier, format="%.2f")
  _, lambda2 = gui.InputDouble("Shrinking ratio", lambda2, format="%.2f")
  changed, outer_radius = gui.DragFloat("Outer radius (mm)", outer_radius, 1, 1, 100, "%.0f")
  if changed and outer_radius > 0:
    P *= 2 * outer_radius / (np.max(P) - np.min(P))
    ps.get_surface_mesh("Parameterization").update_vertex_positions(P)
  changed, eccentricity = gui.DragFloat("Ellipse eccentricity", eccentricity, 1 / 300., 0.01, 1., "%.2f")
  if changed:
    e = (np.max(P[:,1]) - np.min(P[:,1])) / (np.max(P[:,0]) - np.min(P[:,0]))
    P[:,1] *= eccentricity / e
    ps.get_surface_mesh("Parameterization").update_vertex_positions(P)
  _, angle = gui.DragFloat("Spiral pitch angle (degrees)", angle, 1, 0, 90, "%.2f")

  if gui.Button("Generate spiral pattern"):
    # project to disk
    P[:,1] /= eccentricity
    theta = np.arctan2(P[:, 1], P[:, 0])
    theta += np.radians(angle)

    stripe = shrink_morph_py.StripeAlgo(P[:, :2], F)
    layer = stripe.generate_first_layer(P[:, :2], F, theta, printer.nozzle_width)
    nodes, edges = convert_trajectories(layer)

    # rescale
    P[:,1] *= eccentricity
    nodes[:,1] *= eccentricity

    display_trajectories(nodes, edges)

    n_layers = round(thickness / printer.layer_height)
    trajectories = []
    for i in range(n_layers):
      curr_layer = []
      for j in range(len(layer)):
        # rescale
        layer[j][:,1] *= eccentricity

        # rotate wrt previous layer
        s = 2 * printer.nozzle_width / outer_radius
        c = 1 - (2 * printer.nozzle_width / outer_radius)**2 / 2
        rotation_matrix = np.array([[c, -s], [s, c]])
        res = rotation_matrix.dot(layer[j].T)
        layer[j] = res.T

        curr_layer.append(np.column_stack((layer[j], (i + 1) * printer.layer_height * np.ones(layer[j].shape[0]))))
      trajectories.append({"height": printer.layer_height, "paths": curr_layer})

  if gui.Button("Simulate"):
    vertices = []
    segments = []

    make_circle(outer_radius, 1)

    A = dict(vertices=vertices, segments=segments)
    B = tr.triangulate(A, 'pqa1')

    zeros = np.zeros((len(B['vertices']), 1))
    P = np.hstack((np.array(B['vertices']), zeros))
    F = np.array(B['triangles'])

    V = P.copy()
    V[:,2] = 1e-3 * np.random.rand(V.shape[0])

    theta = np.arctan2(P[:, 1] / eccentricity, P[:, 0])
    theta += np.radians(angle)

    shrink_morph_py.simulation(V, P[:,:2], F, theta, 10, lambda1, lambda2, thickness, 1000, 1e-6)

    ps.get_surface_mesh("Parameterization").set_enabled(False)
    ps.register_surface_mesh("Simulation", V, F, edge_width=1, color=(42/255, 53/255, 213/255))


  if gui.Button("Export gcode"):
    filename = filedialog.asksaveasfilename(defaultextension='.gcode')
    printer.to_gcode(trajectories, filename, variable_layer_height=True)
    # printer.to_gcode(trajectories, "output.gcode", variable_layer_height=True)



ps.set_user_callback(callback)


ps.show()