import shrink_morph_py
import numpy as np
import polyscope as ps
import polyscope.imgui as gui
import togcode
import tkinter as tk
from tkinter import filedialog
import os
import sys
import threading

root = tk.Tk()
root.withdraw()

class ShrinkMorph:
  lambda1 = 0.6912916667
  lambda2 = 1.055333333
  gradient = 0
  wD = 2e-5
  E1 = 10
  # deltaLambda = 0.0226764665509417
  n_layers = 8
  lim = 1e-6
  n_iter = 1000
  flattest_print = 0
  width = 200
  wM = 0.01
  wL = 0.01
  thickness = 0.8
  thickness_sim = thickness / lambda1
  num_rectangles = 3
  rect_length = 80
  rect_width = 20
  with_smoothing = False
  printer_profile = "Bambulab_P1S"
  printers_list = [
    'Anet_A8', 
    'Anycubic_Mega_Zero',
    'Artillery_SW_X1', 
    'Bambulab_P1S',
    'Bambulab_X1C',
    'Creality_K1_Max', 
    'Creality_CR-10',
    'Creality_CR-10S_Pro', 
    'Creality_CR-10S5', 
    'E3D_Toolchanger', 
    'Creality_Ender_2', 
    'Creality_Ender_3', 
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
  in_calibration_loop = False
  file_selected = False
  file_not_select_error = False

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

      self.P = shrink_morph_py.parameterization(self.V, self.F, self.lambda1, self.lambda2, self.wD if self.with_smoothing else 0, self.n_iter, self.lim)
      scale = self.width / (np.max(self.P) - np.min(self.P))
      self.P *= scale
      self.V *= scale

      ps.register_surface_mesh("Parameterization", self.P, self.F, material="flat")

      sigma1, sigma2, self.angles = shrink_morph_py.compute_SVD_data(self.V, self.P, self.F)
      ps.get_surface_mesh("Parameterization").add_scalar_quantity("stretch orientation", self.angles, defined_on='faces', enabled=True, vminmax=(-np.pi/2, np.pi/2), cmap='twilight')
      ps.get_surface_mesh("Parameterization").add_scalar_quantity("sigma1", sigma1, defined_on='faces')
      ps.get_surface_mesh("Parameterization").add_scalar_quantity("sigma2", sigma2, defined_on='faces')

      self.file_selected = True
      self.file_not_select_error = False


    changed = gui.BeginCombo("Select printer", self.printer_profile.replace('_', ' '))
    if changed:
      for val in self.printers_list:
        _, selected = gui.Selectable(val.replace('_', ' '), self.printer_profile==val)
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
    if changed and ps.has_surface_mesh("Parameterization"):
      self.P = shrink_morph_py.reparameterization(self.V, self.P, self.F, self.lambda1, self.lambda2, self.wD if self.with_smoothing else 0, self.n_iter, self.lim)
      ps.get_surface_mesh("Parameterization").update_vertex_positions(self.P)

      sigma1, sigma2, self.angles = shrink_morph_py.compute_SVD_data(self.V, self.P, self.F)
      ps.get_surface_mesh("Parameterization").add_scalar_quantity("stretch orientation", self.angles, defined_on='faces', enabled=True, vminmax=(-np.pi/2, np.pi/2), cmap='twilight')
      ps.get_surface_mesh("Parameterization").add_scalar_quantity("sigma1", sigma1, defined_on='faces')
      ps.get_surface_mesh("Parameterization").add_scalar_quantity("sigma2", sigma2, defined_on='faces')

    if ps.has_surface_mesh("Parameterization"):
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
      changed, self.lambda1 = gui.InputDouble("lambda1", self.lambda1, 0, 0, "%.2f")
      changed, self.lambda2 = gui.InputDouble("lambda2", self.lambda2, 0, 0, "%.2f")
      changed, self.thickness_sim = gui.InputDouble("Thickness", self.thickness_sim, 0, 0, "%.2f")
      _, self.gradient = gui.InputDouble("Layer height delta", self.gradient, format="%.3f")
      changed, self.n_iter = gui.InputInt("Iterations", self.n_iter, step=1)
      changed, self.lim = gui.InputDouble("Limit", self.lim, 0, 0, "%.1e")
      gui.TreePop()

    if self.file_not_select_error:
      gui.TextColored((1.0, 0.2, 0.2, 1.0), 
                "(ERROR) No file selected. Please select a file before proceeding.")

    if gui.Button("Next"):
      if not self.file_selected:
        self.file_not_select_error = True
      else:
        self.leave = False
        self.in_calibration_loop = False
        ps.unshow()
      return

    if gui.Button("Calibrate"):
      self.leave = False
      self.calibrate = True
      self.in_calibration_loop = True
      ps.unshow()

    # if gui.Button("Test"):
    #   # V, P, F, theta2 = shrink_morph_py.subdivide(self.V, self.P, self.F, 0 * self.angles, 0.8)
    #   theta1 = shrink_morph_py.vertex_based_stretch_angles(self.V, self.P, self.F)
    #   self.stripe = shrink_morph_py.StripeAlgo(self.P, self.F)
    #   trajectories = self.stripe.generate_one_layer(self.P, self.F, theta1, 0 * theta1, self.printer.layer_height, self.printer.nozzle_width, 10, 0)
    #   # dirs = shrink_morph_py.vertex_based_stretch_vectors(V, P, F)
    #   # trajectories = shrink_morph_py.generate_trajectories(P, F, dirs, 0.4)
    #   nodes, edges = self.convert_trajectories(trajectories)
    #   ps_traj = ps.register_curve_network("Test", nodes, edges, enabled=True, radius=0.2)
    #   ps_traj.set_radius(0.2, relative=False)
    #   ps_traj.set_color((0.3, 0.6, 0.8))
  
  def show(self):
    ps.set_give_focus_on_show(True)
    ps.init()
    ps.set_up_dir("neg_y_up")

    twilight_image = os.path.join(os.path.dirname(__file__), "twilight_colormap.png")
    ps.load_color_map("twilight", twilight_image)
    ps.set_ground_plane_mode("shadow_only")

    # self.V = V
    # self.F = F

    self.param_screen()
    while self.in_calibration_loop:

      if self.leave:
        return
      
      if self.calibrate:
        print("calibrate")
        self.calibrate_screen()

      if self.leave:
        return

      self.after_calibrate_screen()

      if self.leave:
        return

      self.param_screen()

    if self.leave:
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
    changed, self.width = gui.InputDouble("Width", self.width, 0, 0, "%.0f")
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

  curr_z = 0
  
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
    self.n_layers = round(self.thickness / self.printer.layer_height)
    layer = self.stripe.generate_one_layer(self.P, self.F, self.theta1, self.theta2, self.printer.layer_height, self.printer.nozzle_width, self.n_layers, 0)
    nodes, edges = self.convert_trajectories(layer)

    layer_height = self.modified_layer_height(self.printer.layer_height, 0, 1, self.n_layers, self.gradient)
    self.curr_z += layer_height
    for i in range(len(layer)):
      layer[i] = np.column_stack((layer[i], self.curr_z * np.ones(layer[i].shape[0])))
    print(self.curr_z)

    self.trajectories = [{"height": layer_height, "paths": layer}]
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
    spacing = 6
    x_start = -(4 * (self.rect_width+0.1))//2
    step_size = self.rect_width
    build_vert = np.empty(shape=(self.num_rectangles*4, 3))
    build_face = np.empty(shape=(self.num_rectangles, 4))
    y_bottom = -(self.rect_width+spacing) * (self.num_rectangles/2)
    for i in range(self.num_rectangles):
      build_vert[i*4] = [x_start, y_bottom, 0.1]  # Bottom-left
      build_vert[i*4+1] = [x_start, y_bottom+step_size, 0.1]  # Bottom-right
      build_vert[i*4+2] = [x_start+self.rect_length, y_bottom+step_size, 0.1]  # Top-right
      build_vert[i*4+3] = [x_start+self.rect_length, y_bottom, 0.1]  # Top-left
      build_face[i] = [i*4, i*4+1, i*4+2, i*4+3]
      y_bottom += step_size + spacing  # Move upwards for the next rectangle
    ps.register_surface_mesh("Rectangles", build_vert, build_face, color=(0.6, 0.6, 0.3), edge_width=5, edge_color=(0.8, 0.8, 0.8), material="flat")

    ps.set_user_callback(self.callback_calibrate)
    ps.reset_camera_to_home_view()
    ps.show()

  delta = 0.005
  
  def callback_calibrate(self):
    gui.PushItemWidth(200)
    #self.display_buildplate()

    changed = gui.BeginCombo("Select printer", self.printer_profile.replace('_', ' '))
    if changed:
      for val in self.printers_list:
        _, selected = gui.Selectable(val.replace('_', ' '), self.printer_profile==val)
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
    changed_2, self.thickness = gui.InputDouble("Total thickness", self.thickness, format="%.2f")    
    changed_3, self.rect_width = gui.DragFloat("Rectangle width (mm)", self.rect_width, 1, 1, (self.printer.bed_size[1] + 10) / self.num_rectangles - 20, "%.0f")
    changed_4, self.rect_length = gui.DragFloat("Rectangle length (mm)", self.rect_length, 1, 1, self.printer.bed_size[0] - 20, "%.0f")

    if changed_1 or changed_2 or changed_3 or changed_4:
      x_start = -(4 * (self.rect_width+0.1))//2
      step_size = self.rect_width
      spacing = 6
      build_vert = np.empty(shape=(self.num_rectangles*4, 3))
      build_face = np.empty(shape=(self.num_rectangles, 4))
      y_bottom = -(self.rect_width+spacing) * (self.num_rectangles/2)
      for i in range(self.num_rectangles):
        build_vert[i*4] = [x_start, y_bottom, 0.1]  # Bottom-left
        build_vert[i*4+1] = [x_start, y_bottom+step_size, 0.1]  # Bottom-right
        build_vert[i*4+2] = [x_start+self.rect_length, y_bottom+step_size, 0.1]  # Top-right
        build_vert[i*4+3] = [x_start+self.rect_length, y_bottom, 0.1]  # Top-left
        build_face[i] = [i*4, i*4+1, i*4+2, i*4+3]
        y_bottom += step_size + spacing  # Move upwards for the next rectangle
        print(self.rect_length)
        print(self.rect_width)

      ps.register_surface_mesh("Rectangles", build_vert, build_face, color=(0.6, 0.6, 0.3), edge_width=5, edge_color=(0.8, 0.8, 0.8), material="flat")

    _, self.printer.layer_height = gui.InputDouble("Layer Height", self.printer.layer_height, format="%.2f")
    _, self.printer.bed_temp = gui.InputDouble("Bed temperature", self.printer.bed_temp, format="%.0f")
    _, self.printer.extruder_temp = gui.InputDouble("Nozzle temperature", self.printer.extruder_temp, format="%.0f")
    _, self.printer.print_speed = gui.InputDouble("Printing speed (mm/s)", self.printer.print_speed, format="%.0f")
    _, self.printer.first_layer_speed = gui.InputDouble("First layer speed (mm/s)", self.printer.first_layer_speed, format="%.0f")
    _, self.delta = gui.InputDouble("Layer height delta", self.delta, format="%.3f")
    _, self.printer.nloops = gui.InputInt("Skirt loops", self.printer.nloops, step=1)

    if gui.Button("Generate Calibration G-code"):
      layer_height = self.printer.layer_height
      self.generate_calibration(self.num_rectangles, self.printer.layer_height, self.thickness, self.printer.nozzle_width, self.rect_length, self.rect_width, self.delta)
      self.printer.layer_height = layer_height
    if gui.Button("Finished Calibration"):
      self.leave = False
      ps.unshow()

    # if gui.Button("Back"):
    #     self.calibrate = False
    #     ps.unshow()
    #     self.param_screen()
  def after_calibrate_screen(self):
    self.leave = True
    ps.remove_all_structures()

    self.display_buildplate()
    spacing = 6
    ps.set_user_callback(self.callback_after_calibrate)
    ps.reset_camera_to_home_view()
    ps.show()

  delta = 0.005
  
  def callback_after_calibrate(self):
    gui.PushItemWidth(200)
    #self.display_buildplate()

    changed = gui.BeginCombo("Select printer", self.printer_profile.replace('_', ' '))
    if changed:
      for val in self.printers_list:
        _, selected = gui.Selectable(val.replace('_', ' '), self.printer_profile==val)
        if selected:
          self.printer_profile = val
          self.printer = togcode.Printer(self.printer_profile)
          build_vert = np.array([[-1,-1,-0.1], [1, -1, -0.1], [1, 1, -0.1], [-1, 1, -0.1]])
          build_vert[:, 0] *= self.printer.bed_size[0] / 2
          build_vert[:, 1] *= self.printer.bed_size[1] / 2
          ps.get_surface_mesh("Buildplate").update_vertex_positions(build_vert)
      gui.EndCombo()
    gui.PopItemWidth()

    x_start = -(4 * (self.rect_width+0.1))//2
    step_size = self.rect_width
    spacing = 6
    build_vert = np.empty(shape=(self.num_rectangles*4, 3))
    build_face = np.empty(shape=(self.num_rectangles, 4))
    y_bottom = -(self.rect_width+spacing) * (self.num_rectangles/2)
    for i in range(self.num_rectangles):
      build_vert[i*4] = [x_start, y_bottom, 0.1]  # Bottom-left
      build_vert[i*4+1] = [x_start, y_bottom+step_size, 0.1]  # Bottom-right
      build_vert[i*4+2] = [x_start+self.rect_length, y_bottom+step_size, 0.1]  # Top-right
      build_vert[i*4+3] = [x_start+self.rect_length, y_bottom, 0.1]  # Top-left
      build_face[i] = [i*4, i*4+1, i*4+2, i*4+3]
      y_bottom += step_size + spacing  # Move upwards for the next rectangle
      print(self.rect_length)
      print(self.rect_width)

    ps.register_surface_mesh("Rectangles", build_vert, build_face, color=(0.6, 0.6, 0.3), edge_width=5, edge_color=(0.8, 0.8, 0.8), material="flat")
    changed, self.flattest_print = gui.InputInt("Select the flattest print", self.flattest_print, step=1)
    
    if gui.Button("Confirm"):
      self.leave = False
      ps.reset_camera_to_home_view()
      ps.remove_all_structures()
      ps.unshow()

    gui.PushItemWidth(100)

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
      layer = self.stripe.generate_one_layer(self.P, self.F, self.theta1, self.theta2, self.printer.layer_height, self.printer.nozzle_width, self.n_layers, self.curr_layer)
      
      layer_height = self.modified_layer_height(self.printer.layer_height, self.curr_layer, 1, self.n_layers, self.gradient)
      self.curr_z += layer_height
      for i in range(len(layer)):
        layer[i] = np.column_stack((layer[i], self.curr_z * np.ones(layer[i].shape[0])))
      print(self.curr_z)

      self.trajectories.append({"height": layer_height, "paths": layer})
      nodes, edges = self.convert_trajectories(layer)
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
      self.printer.to_gcode(self.trajectories, filename, variable_layer_height=True)
  
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
  
  def modified_layer_height(self, layer_height, layer_id, rectangle_id, nb_layers, gradient):
    return layer_height + rectangle_id * (layer_id / (nb_layers - 1) - 1 / 2) * gradient

  def generate_calibration(self, nb_rectangles, layer_height, total_thickness, nozzle_width, rect_length, rect_width, gradient, gap_width=10):
    nb_layers = round(total_thickness / layer_height)
    layers = []
    posZ = np.zeros(nb_rectangles)
    for i in range(nb_layers):
        posY = 0
        for j in range(nb_rectangles):
          posY -= rect_width + gap_width
          height = self.modified_layer_height(layer_height, i, j - (nb_rectangles - 1) / 2, nb_layers, gradient)
          posZ[j] += height
          layer = {"height": height, 
                   "paths": self.zigzag_layer(rect_length - self.printer.nozzle_width, rect_width, posY, posZ[j], nozzle_width)}
          layers.append(layer)
          print(f"layer height: {height:.4f}; posZ: {posZ[j]:.4f}") # for debug purposes
    save_path = filedialog.asksaveasfilename(defaultextension="gcode", initialdir=os.getcwd())
    self.printer.to_gcode(layers, save_path, variable_layer_height=True)

main = ShrinkMorph()
main.show()