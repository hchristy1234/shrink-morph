import shrink_morph_py
import numpy as np
import polyscope as ps
import polyscope.imgui as gui
import togcode
import tkinter as tk
from tkinter import filedialog
import os
import sys

root = tk.Tk()
root.withdraw()

class ShrinkMorph:
  lambda1 = 0.58
  lambda2 = 1.08
  wD = 1e-5
  E1 = 10
  thickness_sim = 1.218
  # deltaLambda = 0.0226764665509417
  n_layers = 10
  lim = 1e-6
  n_iter = 1000
  width = 200
  wM = 0.01
  wL = 0.01
  thickness = 0.8
  num_rectangles = 5
  num_layers = 10
  layer_height = 0.08
  rect_length = 80
  rect_width = 20
  with_smoothing = False
  printer_profile = "Prusa_MK3S"
  printers_list = [
    'Anet_A8', 
    'Anycubic_Mega_Zero',
    'Artillery_SW_X1', 
    'BambuLab_X1C',
    'Creality_K1_Max', 
    'CR10',
    'CR10S_Pro', 
    'CR10_S5', 
    'E3D_Toolchanger', 
    'Ender2', 
    'Ender3', 
    'Flsun_SR', 
    'Kingroon_P3', 
    'Kywoo3D_Tycoon'
    'Lutum_4.6', 
    'Micro_Delta_Rework', 
    'Monoprice_select_mini_v2', 
    'Pharaoh_xd', 
    'Prusa_MK2S', 
    'Prusa_MK3S', 
    'Replicator2', 
    'Replicator2x', 
    'Snapmaker', 
    'Strateo3D', 
    'Ultimaker2', 
    'Ultimaker3', 
    'UltimakerS3', 
    'UltimakerS5', 
    'UltimakerS7', 
    'Voron_V0', 
    'Volumic_Stream30_Pro_MK2', 
    'Wasp_2040_Pro', 
  ]
  printer = togcode.Printer(printer_profile)

  # GUI variables
  leave = True
  calibrate = False

  # Display printer buildplate
  def display_buildplate(self):
    build_vert = np.array([[-1,-1,-0.1], [1, -1, -0.1], [1, 1, -0.1], [-1, 1, -0.1]])
    build_vert[:, 0] *= self.printer.bed_size[0] / 2
    build_vert[:, 1] *= self.printer.bed_size[1] / 2
    build_face = np.array([[0, 1, 2, 3]])

    ps.register_surface_mesh("Buildplate", build_vert, build_face, color=(0.95, 0.95, 0.95), edge_width=5, edge_color=(0.5, 0.5, 0.5), material="flat")

  def param_screen(self):
    print("in param screen")

    self.display_buildplate()

    ps.set_user_callback(self.callback_param)
    ps.show()

  def callback_param(self):
    # global lambda1, lambda2, wD, E1, thickness_sim, deltaLambda, n_layers, lim, n_iter, width, P, self.V, self.F, printer_profile, with_smoothing
    gui.PushItemWidth(200)

    if gui.Button("Select File"):
      filename = filedialog.askopenfilename(defaultextension=".obj", filetypes=[("Wavefront OBJ", "*.obj")])
      self.V, self.F = shrink_morph_py.read_from_OBJ(filename)

      self.P = shrink_morph_py.parameterization(self.V, self.F, self.lambda1, self.lambda2, 0, self.n_iter, self.lim)
      scale = self.width / (np.max(self.P) - np.min(self.P))
      self.P *= scale
      self.V *= scale

      ps.register_surface_mesh("Parameterization", self.P, self.F, material="flat")

      sigma1, sigma2, self.angles = shrink_morph_py.compute_SVD_data(self.V, self.P, self.F)
      ps.get_surface_mesh("Parameterization").add_scalar_quantity("stretch orientation", self.angles, defined_on='faces', enabled=True, vminmax=(-np.pi/2, np.pi/2), cmap='twilight')
      ps.get_surface_mesh("Parameterization").add_scalar_quantity("sigma1", sigma1, defined_on='faces')
      ps.get_surface_mesh("Parameterization").add_scalar_quantity("sigma2", sigma2, defined_on='faces')


    changed = gui.BeginCombo("Select printer", self.printer_profile)
    if changed:
      for val in self.printers_list:
        _, selected = gui.Selectable(val, self.printer_profile==val)
        if selected:
          self.printer_profile = val
          self.printer = togcode.Printer(self.printer_profile)
          build_vert = np.array([[-1,-1,-0.1], [1, -1, -0.1], [1, 1, -0.1], [-1, 1, -0.1]])
          build_vert[:, 0] *= self.printer.bed_size[0] / 2
          build_vert[:, 1] *= self.printer.bed_size[1] / 2
          ps.get_surface_mesh("Buildplate").update_vertex_positions(build_vert)
      gui.EndCombo()
    gui.PopItemWidth()

    gui.PushItemWidth(100)
    changed, self.with_smoothing = gui.Checkbox("With smoothing", self.with_smoothing) 
    if changed:
      self.P = shrink_morph_py.reparameterization(self.V, self.P, self.F, self.lambda1, self.lambda2, self.wD if self.with_smoothing else 0, self.n_iter, self.lim)
      ps.get_surface_mesh("Parameterization").update_vertex_positions(self.P)

      sigma1, sigma2, self.angles = shrink_morph_py.compute_SVD_data(self.V, self.P, self.F)
      ps.get_surface_mesh("Parameterization").add_scalar_quantity("stretch orientation", self.angles, defined_on='faces', enabled=True, vminmax=(-np.pi/2, np.pi/2), cmap='twilight')
      ps.get_surface_mesh("Parameterization").add_scalar_quantity("sigma1", sigma1, defined_on='faces')
      ps.get_surface_mesh("Parameterization").add_scalar_quantity("sigma2", sigma2, defined_on='faces')
    changed, self.width = gui.DragFloat("Width", self.width, 1, 1, 500, "%.0f")
    if changed and self.width > 0:
      scale = self.width / (np.max(self.P) - np.min(self.P))
      self.P *= scale
      self.V *= scale
      ps.get_surface_mesh("Parameterization").update_vertex_positions(self.P)

    if gui.Button("Increase mesh resolution"):
      self.V, self.P, self.F, _ = shrink_morph_py.subdivide(self.V, self.P, self.F, np.array([]))
      self.P = shrink_morph_py.reparameterization(self.V, self.P, self.F, self.lambda1, self.lambda2, self.wD if self.with_smoothing else 0, self.n_iter, self.lim)
      ps.register_surface_mesh("Parameterization", self.P, self.F, material="flat")

      sigma1, sigma2, self.angles = shrink_morph_py.compute_SVD_data(self.V, self.P, self.F)
      ps.get_surface_mesh("Parameterization").add_scalar_quantity("stretch orientation", self.angles, defined_on='faces', enabled=True, vminmax=(-np.pi/2, np.pi/2), cmap='twilight')
      ps.get_surface_mesh("Parameterization").add_scalar_quantity("sigma1", sigma1, defined_on='faces')
      ps.get_surface_mesh("Parameterization").add_scalar_quantity("sigma2", sigma2, defined_on='faces')

    if gui.TreeNode("Advanced"):
      changed, self.lambda1 = gui.InputFloat("self.lambda1", self.lambda1, 0, 0, "%.1e")
      changed, self.lambda2 = gui.InputFloat("self.lambda2", self.lambda2, 0, 0, "%.1e")
      changed, self.thickness_sim = gui.InputFloat("self.Thickness", self.thickness_sim, 0, 0, "%.1e")
      changed, self.n_iter = gui.InputInt("Iterations", self.n_iter, step=1)
      changed, self.lim = gui.InputFloat("self.Limit", self.lim, 0, 0, "%.1e")
      gui.TreePop()

    if gui.Button("Next"):
      self.leave = False
      ps.unshow()

    if gui.Button("Calibrate"):
      self.leave = False
      self.calibrate = True
      ps.unshow()
  
  def show(self):
    ps.set_give_focus_on_show(True)
    ps.init()
    ps.set_up_dir("neg_y_up")

    twilight_image = os.path.join(os.path.dirname(__file__), "twilight_colormap.png")
    print(twilight_image)
    ps.load_color_map("twilight", twilight_image)

    #ps.load_color_map("twilight", "data/twilight_colormap.png")
    ps.set_ground_plane_mode("shadow_only")

    # self.V = V
    # self.F = F

    self.param_screen()

    if self.leave:
      return
    
    if self.calibrate:
      print("calibrate")
      self.calibrate_screen()
      return

    self.optim_screen()

    if self.leave:
      return
    
    self.traj_screen()


  # Directions optimization
  def optim_screen(self):
    self.leave = True
    ps.remove_all_structures()
    ps_input = ps.register_surface_mesh("Input mesh", self.V, self.F)
    ps_input.set_transparency(0.5)

    self.targetV = self.V.copy()
    self.theta2 = np.zeros(self.V.shape[0])
    self.optim_solver = shrink_morph_py.SGNSolver(self.targetV, self.P, self.F, self.E1, self.lambda1, self.lambda2, self.thickness_sim)

    ps.set_user_callback(self.callback_optim)
    ps.reset_camera_to_home_view()
    ps.show()

  resolutions = ["Low", "Medium", "High"]
  resolution = resolutions[0]

  optim_running = False

  def callback_optim(self):
    gui.PushItemWidth(100)
    changed, self.width = gui.InputFloat("Width", self.width, 0, 0, "%.0f")
    if changed and self.width > 0:
      scale = self.width / (np.max(self.P) - np.min(self.P))
      self.P *= scale
      self.V *= scale
      self.targetV *= scale
      ps.get_surface_mesh("Input mesh").update_vertex_positions(self.targetV)

    # if gui.Button("Simulation"):
    #   shrink_morph_py.simulation(self.V, self.P, self.F, self.theta2, self.E1, self.lambda1, self.lambda2, self.deltaLambda, self.thickness_sim, self.width, self.n_iter, self.lim)
    #   ps.get_surface_mesh("Input mesh").set_transparency(0.5)
    #   ps.register_surface_mesh("Simulation", self.V, self.F)

    if self.optim_running == True:
      _, self.theta2 = self.optim_solver.solve_one_step()
      ps.get_surface_mesh("Optimized mesh").update_vertex_positions(self.optim_solver.optimizedV())

      if self.optim_solver.decrement() < 1e-6:
        self.optim_running = False
        ps.get_surface_mesh("Optimized mesh").add_scalar_quantity("theta2", self.theta2)
        ps.get_surface_mesh("Optimized mesh").add_scalar_quantity("theta1", self.angles, defined_on='faces', vminmax=(-np.pi/2, np.pi/2), cmap='twilight')

    if gui.Button("Directions optimization"):
      ps.get_surface_mesh("Input mesh").set_transparency(0.5)
      # if ps.has_surface_mesh("Simulation"):
      #   ps.get_surface_mesh("Simulation").update_vertex_positions(self.V)
      # else:
      print("Initial distance", self.optim_solver.distance(self.theta2))
      ps.register_surface_mesh("Optimized mesh", self.optim_solver.optimizedV(), self.F)
      self.optim_running = True

    changed = gui.BeginCombo("Trajectory resolution", self.resolution)
    if changed:
      for val in self.resolutions:
        _, selected = gui.Selectable(val, self.resolution==val)
        if selected:
          self.resolution = val
      gui.EndCombo()

    if gui.Button("Generate trajectories"):
      self.leave = False
      ps.unshow()
  
  def traj_screen(self):
    # Trajectories & G-code generation
    if self.resolution == "Low":
      target_edge_length = 1
    elif self.resolution == "Medium":
      target_edge_length = 0.5
    elif self.resolution == "High":
      target_edge_length = 0.2
    self.V, self.P, self.F, self.theta2 = shrink_morph_py.subdivide(self.V, self.P, self.F, self.theta2, target_edge_length)
    self.theta1 = shrink_morph_py.vertex_based_stretch_angles(self.V, self.P, self.F)
    self.stripe = shrink_morph_py.StripeAlgo(self.P, self.F)
    trajectories = self.stripe.generate_one_layer(self.P, self.F, self.theta1, self.theta2, self.printer.layer_height, self.printer.nozzle_width, self.n_layers, 0)
    nodes, edges = self.convert_trajectories(trajectories)
    self.layer_nodes = [nodes]
    self.layer_edges = [edges]
    ps.remove_all_structures()
    self.display_trajectories(self.layer_nodes, self.layer_edges)
  
    ps.set_user_callback(self.callback_traj)
    ps.reset_camera_to_home_view()
    ps.show()
  
  def calibrate_screen(self):
    self.leave = True
    ps.remove_all_structures()

    self.display_buildplate()
    x_start = -(self.num_rectangles * (self.rect_width+0.1))//2
    step_size = self.rect_width
    build_vert = np.empty(shape=(int(self.num_rectangles*4), 3))
    build_face = np.empty(shape=(int(self.num_rectangles), 4))
    y_bottom = -self.rect_length / 2
    y_top = self.rect_length / 2
    for i in range(int(self.num_rectangles)):
      build_vert[i*4] = [x_start, y_bottom, 0.1]
      build_vert[i*4+1] = [x_start+step_size, y_bottom, 0.1]
      build_vert[i*4+2] = [x_start+step_size, y_top, 0.1]
      build_vert[i*4+3] =  [x_start, y_top, 0.1]
      build_face[i] = [i*4, i*4+1, i*4+2, i*4+3]
      x_start += step_size + 0.05
    ps.register_surface_mesh("Rectangles", build_vert, build_face, color=(0.6, 0.6, 0.3), edge_width=5, edge_color=(0.8, 0.8, 0.8), material="flat")

    ps.set_user_callback(self.callback_calibrate)
    ps.reset_camera_to_home_view()
    ps.show()

  gradient = 0.0005
  discrete = False
  
  def callback_calibrate(self):
    gui.PushItemWidth(200)
    #self.display_buildplate()

    changed = gui.BeginCombo("Select printer", self.printer_profile)
    if changed:
      for val in self.printers_list:
        _, selected = gui.Selectable(val, self.printer_profile==val)
        if selected:
          self.printer_profile = val
          self.printer = togcode.Printer(self.printer_profile)
          build_vert = np.array([[-1,-1,-0.1], [1, -1, -0.1], [1, 1, -0.1], [-1, 1, -0.1]])
          build_vert[:, 0] *= self.printer.bed_size[0] / 2
          build_vert[:, 1] *= self.printer.bed_size[1] / 2
          ps.get_surface_mesh("Buildplate").update_vertex_positions(build_vert)
      gui.EndCombo()
    gui.PopItemWidth()

    gui.PushItemWidth(100)

    changed_1, self.num_rectangles = gui.InputInt("Number of Rectangles", self.num_rectangles, step=1)
    changed_2, self.thickness = gui.InputFloat("Total thickness", self.thickness, format="%.2f")    
    changed_3, self.rect_width = gui.DragFloat("Rectangle width (mm)", self.rect_width, 1, 1, (self.printer.bed_size[1] + 10) / self.num_rectangles - 20, "%.0f")
    changed_4, self.rect_length = gui.DragFloat("Rectangle length (mm)", self.rect_length, 1, 1, self.printer.bed_size[0] - 20, "%.0f")

    if changed_1 or changed_2 or changed_3 or changed_4:
      x_start = -(self.num_rectangles * (self.rect_width+0.1))//2
      step_size = self.rect_width
      build_vert = np.empty(shape=(int(self.num_rectangles*4), 3))
      build_face = np.empty(shape=(int(self.num_rectangles), 4))
      y_bottom = -self.rect_length / 2
      y_top = self.rect_length / 2
      for i in range(int(self.num_rectangles)):
        build_vert[i*4] = [x_start, y_bottom, 0.1]
        build_vert[i*4+1] = [x_start+step_size, y_bottom, 0.1]
        build_vert[i*4+2] = [x_start+step_size, y_top, 0.1]
        build_vert[i*4+3] =  [x_start, y_top, 0.1]
        build_face[i] = [i*4, i*4+1, i*4+2, i*4+3]
        x_start += step_size + 0.1

      ps.register_surface_mesh("Rectangles", build_vert, build_face, color=(0.6, 0.6, 0.3), edge_width=5, edge_color=(0.8, 0.8, 0.8), material="flat")

    _, self.printer.layer_height = gui.InputFloat("Layer Height", self.printer.layer_height, format="%.2f")
    _, self.printer.bed_temp = gui.InputFloat("Bed temperature", self.printer.bed_temp, format="%.0f")
    _, self.printer.extruder_temp = gui.InputFloat("Nozzle temperature", self.printer.extruder_temp, format="%.0f")
    _, self.printer.print_speed = gui.InputFloat("Printing speed (mm/s)", self.printer.print_speed, format="%.0f")
    _, self.printer.first_layer_speed = gui.InputFloat("First layer speed (mm/s)", self.printer.first_layer_speed, format="%.0f")
    _, self.gradient = gui.InputFloat("Gradient", self.gradient, format="%.4f")
    _, self.discrete = gui.Checkbox("Discrete", self.discrete)
    _, self.printer.nloops = gui.InputInt("Skirt loops", self.printer.nloops, step=1)

    if gui.Button("Generate Calibration G-code"):
      layer_height = self.printer.layer_height
      self.generate_calibration(self.num_rectangles, self.printer.layer_height, self.thickness, self.printer.nozzle_width, self.rect_length, self.rect_width, self.gradient, self.discrete)
      self.printer.layer_height = layer_height
    # if gui.Button("Back"):
    #     self.calibrate = False
    #     ps.unshow()
    #     self.param_screen()

  def convert_trajectories(self, trajectories):
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

  def display_trajectories(self, nodes, edges):
    self.display_buildplate()
    for k in range(len(nodes)):
      ps_traj = ps.register_curve_network("Layer " + str(k + 1), nodes[k], edges[k], enabled=True, radius=self.printer.nozzle_width / 2)
      ps_traj.set_radius(self.printer.nozzle_width / 2, relative=False)
      ps_traj.set_color((0.3, 0.6, 0.8))

  layer_id = 1
  progress = 1
  curr_layer = 1
  def callback_traj(self):
    if self.curr_layer < self.n_layers:
      trajectories = self.stripe.generate_one_layer(self.P, self.F, self.theta1, self.theta2, self.printer.layer_height, self.printer.nozzle_width, self.n_layers, self.curr_layer)
      nodes, edges = self.convert_trajectories(trajectories)
      self.layer_nodes.append(nodes)
      self.layer_edges.append(edges)
      self.display_trajectories(self.layer_nodes, self.layer_edges)
      self.curr_layer = self.curr_layer + 1

    gui.PushItemWidth(200)
    changed = gui.BeginCombo("Select printer", self.printer_profile)
    if changed:
      for val in self.printers_list:
        _, selected = gui.Selectable(val, self.printer_profile==val)
        if selected:
          self.printer_profile = val
          self.printer = togcode.Printer(self.printer_profile)
          self.display_buildplate()
      gui.EndCombo()
    gui.PopItemWidth()


    gui.PushItemWidth(200)
    changed, self.layer_id = gui.SliderInt("Layer", self.layer_id, 1, self.curr_layer)
    if changed:
      for i in range(1, self.n_layers + 1):
        if i == self.layer_id:
          ps.get_curve_network("Layer " + str(i)).set_enabled(True)
        else:
          ps.get_curve_network("Layer " + str(i)).set_enabled(False)
    changed, self.progress = gui.SliderInt("Progress", self.progress, 1, self.layer_edges[self.layer_id - 1].shape[0])
    gui.PopItemWidth()
    if changed:
      d = int(np.max(self.layer_edges[self.layer_id - 1][:self.progress, :])) + 1
      ps.register_curve_network("Layer " + str(self.layer_id), self.layer_nodes[self.layer_id - 1][:d, :], self.layer_edges[self.layer_id - 1][:self.progress, :])
    if gui.Button("Increase mesh resolution and reload trajectories"):
      self.V, self.P, self.F, self.theta2 = shrink_morph_py.subdivide(self.V, self.P, self.F, self.theta2)
      self.theta1 = shrink_morph_py.vertex_based_stretch_angles(self.V, self.P, self.F)
      self.stripe = shrink_morph_py.StripeAlgo(self.P, self.F)
      self.layer_nodes = []
      self.layer_edges = []
      self.curr_layer = 0
      ps.remove_all_structures()
    if gui.Button("Export to g-code"):
      filename = filedialog.asksaveasfilename(defaultextension='.gcode')
      self.printer.to_gcode(self.trajectories, filename)
  
  def read_trajectories(self, filename):
    print("reading file " + filename)

    with open(filename, "r") as file:
        paths = []
        line = file.readline()
        while line:
            num_vertices = int(line)
            path = []
            for _ in range(num_vertices):
                line = file.readline()
                vertex = line.split()
                path.append([float(x) for x in vertex])
            paths.append(np.array(path))
            line = file.readline()
    return paths

  def zigzag_layer(self, length, width, posY, posZ, nozzle_width):
    left = True
    y = width / 2
    paths = []
    while(y > -width / 2):
      if left: # go from left to right
        paths.append(np.array([[0, posY + y, posZ], [length, posY + y, posZ]]))
        left = False
      else: # go from right to left
        paths.append(np.array([[length, posY + y, posZ], [0, posY + y, posZ]]))
        left = True
      y -= nozzle_width
    return paths
  
  def modified_layer_height(self, layer_height, layer_id, rectangle_id, nb_layers, gradient, discrete):
      if discrete:
        if layer_id == 0:
          return layer_height - rectangle_id * gradient
        else:
          return layer_height + rectangle_id * gradient / (nb_layers - 1)
      else:
        return layer_height + rectangle_id * (layer_id - (nb_layers - 1) / 2) * gradient

  def generate_calibration(self, nb_rectangles, layer_height, total_thickness, nozzle_width, rect_length, rect_width, gradient, discrete, gap_width=10):
    nb_layers = round(total_thickness / layer_height)
    paths = []
    posZ = np.zeros(nb_rectangles)
    for i in range(nb_layers):
        posY = 0
        for j in range(nb_rectangles):
          posY -= rect_width + gap_width
          posZ[j] += self.modified_layer_height(layer_height, i, j, nb_layers, gradient, discrete)
          paths = paths + self.zigzag_layer(rect_length - self.printer.nozzle_width, rect_width, posY, posZ[j], nozzle_width)
          print(f"{posZ[j]:.4f}") # for debug purposes
    save_path = filedialog.asksaveasfilename(defaultextension="gcode", initialdir=os.getcwd())
    self.printer.to_gcode(paths, save_path)

main = ShrinkMorph()
main.show()